# -*- coding: utf-8 -*-
import time
import os,sys,glob,copy,pickle
import torch
import torchvision
import torchvision.transforms as transforms
import random
from build_fg_cpu import build_data_strict_all,build_gkernel,\
                     build_data_smtry,build_data_strict_smtry,build_data_strict_paths
import numpy as np
import torch #new 230518
import torch.nn as nn
import torch.nn.functional as F
import networks
from scipy import ndimage #new
import torch.optim as optim

def write_out_pdb(fpath,vecs,scores,scorecut = None):
    form = 'HETATM%5d  O   HOH X%4d    %8.3f%8.3f%8.3f  1.00%6.2f\n'
    length = len(vecs)
    f = open(fpath,'w')
    dump = []
    for i in range(length):
        write = False
        if (scorecut == None):
            write = True
        elif scores[i] > scorecut:
            write = True
        else:
            write = False
        if write:
            newline = form%((i+1),(1),vecs[i][0],vecs[i][1],vecs[i][2],scores[i])
            dump.append( [ scores[i] , newline ] )
    #sort!
    dump_sorted = sorted(dump, key = lambda x: -1.0*x[0])
    for i in range(len(dump_sorted)):
        f.write(dump_sorted[i][1])
    f.close()
def place_water(prob, null, vec_start, kernel, vcut=1.0, grid =0.5,padding=4.0,ignore_padding=False):
    #prob, ban: 3 dimension array
    #probability: result from net
    #null: True when atom(C,N,O, water) not exists 
    #vec_start: coord of prob[0][0][0]

    #build oxygen sized gaussian kernel
    #1.52: radius of oxygen atom (in angstrom) 
    r_grid = int( 1 + (1.5*1.52)/grid ) #radius in grid debug
    #r_grid = int( 1 + (1.52)/grid ) #radius in grid
    d_grid = 2*r_grid+1
    result = [] 
    scores = [] 
    maxv = 999.999
   
    #new 
    prob_new  = ndimage.convolve(prob, kernel, mode='reflect')
    cv = prob_new * null
    #faster algorithm
    # 1. sort for cv_flatten, which is flattened version of origial array, cv. 
    # 2. use argsort for finding index of highest scored water
    # 3. check cv[idx] > vcut. 
    # 4. place water
    # 5. fill negative value to the too-close(<2.28A) element of cv array.
    device = torch.device("cuda:0" if torch.cuda.is_available() else "cpu")
    vec_init = vec_start #vec_start for ignore_padding is False, vec_start+padding for ignore_padding is True.
    #removeing padded area
    if ignore_padding:
        vec_init += np.array([padding,padding,padding])
        eps = 0.0001#for numerical stability        
        pad_idx = int((padding +eps) / grid)
        cv = cv[pad_idx:-pad_idx,pad_idx:-pad_idx,pad_idx:-pad_idx]
   
    cv_torch = torch.tensor(cv, dtype=torch.float32, device=device, requires_grad=False)
    cv_flatten = torch.flatten(cv_torch)
    cv_argsort = torch.argsort(cv_flatten,descending=True) 
    cv_argsort_np = cv_argsort.cpu().detach().numpy()
    cv_shape = cv.shape
    n_cv = cv_flatten.shape[0]
    for wat_idx in range(n_cv):
        #finding index
        argmax    = cv_argsort[wat_idx]
        argmax_np = cv_argsort_np[wat_idx]
        # using // for torch has further compatibility issue
        ind    = ( torch.div(argmax, (cv_shape[1]*cv_shape[2]), rounding_mode='floor'),
                   torch.div(argmax,  cv_shape[2]             , rounding_mode='floor')%cv_shape[1],
                   argmax%cv_shape[2] )        
        ind_np = ( argmax_np//(cv_shape[1]*cv_shape[2]), (argmax_np//cv_shape[2])%cv_shape[1], argmax_np%cv_shape[2] )
        val = cv_torch[ind]

        if val <0: #too-close from placed water
            continue
        elif val < vcut:
            break  

        vec = vec_init + grid*np.array(ind_np)
        grid_vec = grid*torch.tensor(ind,dtype=torch.float32, device=device, requires_grad=False)#for water removal purpose
        result.append(vec)
        scores.append(val.cpu().detach().numpy())
        #remove water score from CV too close with water elements
        #check btorch.py
        
        min_pt = [ max(0, ind[i]-r_grid) for i in range(3)]
        max_pt = [ min(cv.shape[i]-1, ind[i]+r_grid) for i in range(3)]        
        x = torch.arange(min_pt[0],max_pt[0],1, out=torch.cuda.LongTensor() )
        y = torch.arange(min_pt[1],max_pt[1],1, out=torch.cuda.LongTensor() )
        z = torch.arange(min_pt[2],max_pt[2],1, out=torch.cuda.LongTensor() )
        xind,yind,zind =torch.meshgrid([x,y,z],indexing='ij') #add indexing
        xv = grid*xind.float()
        yv = grid*yind.float()
        zv = grid*zind.float()
        vec_x = torch.full_like(xv,grid_vec[0])
        vec_y = torch.full_like(yv,grid_vec[1])
        vec_z = torch.full_like(zv,grid_vec[2])
        #
        dx = vec_x-xv
        dy = vec_y-yv
        dz = vec_z-zv
        #
        dx2 = dx*dx 
        dy2 = dy*dy 
        dz2 = dz*dz 
        #
        d2 = (dx2+dy2+dz2)
        d  = torch.sqrt(d2)
        #torch.where: (contidion, true,false)
        mask  = torch.where(d<(1.5*1.52),-1000.0,0.0) #arbitary value big enough to make cv negative
        cv_torch[min_pt[0]:max_pt[0], min_pt[1]:max_pt[1], min_pt[2]:max_pt[2] ] += mask[:,:,:]        
    return result,scores

def batch_idxs(idxs, batch=4, shuffle=True):
    result = []
    if shuffle:
        random.shuffle(idxs)
    len_idxs = len(idxs)
    for i in range(len_idxs//batch):
        result.append([idxs[batch*i+j] for j in range(batch)])
    return result

def build_batch(batch,train=True,n_grid=32,dbloss=False):
    inputs_l = []
    outputs_l = []
    w_l = []
    vec_start = []
    for idx in batch:
        if train:
            vs,inp,out,w = build_data_smtry(idx,n_grid=n_grid,dbloss=dbloss)
        else:
            vs,inp,out,w = build_data_strict_smtry(idx,n_grid=n_grid,dbloss=dbloss)
        inputs_l.append(inp)
        outputs_l.append(out)
        w_l.append(w)
        vec_start.append(vs)
    
    inputs  = np.vstack(inputs_l)
    outputs = np.vstack(outputs_l)
    ws = np.vstack(w_l)
    return vec_start,inputs,outputs,ws       

def build_batch_all(batch,padding=32.0,dbloss=False):
    inputs_l = []
    outputs_l = []
    w_l = []
    vec_start = []
    for idx in batch:
        vs,inp,out,w = build_data_strict_all(idx,padding=padding,dbloss=dbloss)
        inputs_l.append(inp)
        outputs_l.append(out)
        w_l.append(w)
        vec_start.append(vs)
    inputs  = np.vstack(inputs_l)
    outputs = np.vstack(outputs_l)
    ws = np.vstack(w_l)
    return vec_start,inputs,outputs,ws       

def build_batch_full_paths(paths,grid=0.5,n_grid=64,padding=4.0,dbloss=False):
    #'paths : {"pro":[pro_paths...],"wat":[wat_paths(=answer,for training, optional)] }
    inputs_l = []
    outputs_l = []
    w_l = []
    vec_start = []
    data_old = build_data_strict_paths(paths,grid=0.5,n_grid=64,padding=4.0,dbloss=dbloss)
    data = []
    for datum in data_old:
        vs  = np.array( [datum[0]] )
        inp = np.array( datum[1] )
        out = np.array( datum[2] )
        w   = np.array( datum[3] )
        data.append( (vs,inp,out,w))
    return data       

def save_grid(vec_starts, np_inputs, out, batch, dr ,vcut=0.6,grid=0.5 ):
    kernel = build_gkernel(grid=grid, r = 1.52)
    np_o    = np.expand_dims((np_inputs[:,0,:,:,:] > 0.5),axis=1)  
    np_c    = np.expand_dims((np_inputs[:,1,:,:,:] > 0.5),axis=1) 
    np_n    = np.expand_dims((np_inputs[:,2,:,:,:] > 0.5),axis=1) 
    np_s    = np.expand_dims((np_inputs[:,3,:,:,:] > 0.5),axis=1) 
    np_prot = np.logical_or(np.logical_or(np.logical_or(np_o,np_c),np_n),np_s)
    np_null = np.logical_not(np_prot)
    
    for trg_idx in range((out.shape[0])): 
        zipped = {'vec_start':vec_starts[trg_idx],
                  'wat_grid':out[trg_idx][0],
                  'wat_grid_noprot':(out[trg_idx][0] * np_null[trg_idx][0]),
                  'grid':grid}
        
        fpath = '%s/%s_grid.bin'%(dr,batch[trg_idx])
        grid_dl = open(fpath,'wb')
        pickle.dump(zipped,grid_dl)
def run_place_water(vec_starts, np_inputs, out, batch, dr ,vcut=0.6,grid=0.5 ):
    kernel = build_gkernel(grid=grid, r = 1.52) 
    np_o    = np.expand_dims((np_inputs[:,0,:,:,:] > 0.5),axis=1)  
    np_c    = np.expand_dims((np_inputs[:,1,:,:,:] > 0.5),axis=1) 
    np_n    = np.expand_dims((np_inputs[:,2,:,:,:] > 0.5),axis=1) 
    np_s    = np.expand_dims((np_inputs[:,3,:,:,:] > 0.5),axis=1) 
    np_prot = np.logical_or(np.logical_or(np.logical_or(np_o,np_c),np_n),np_s)
    np_null = np.logical_not(np_prot)
    
    for trg_idx in range((out.shape[0])): 
        watvecs,scores =  place_water(out[trg_idx][0], np_null[trg_idx][0], 
                                    vec_starts[trg_idx], kernel, vcut=vcut, grid =grid)
        zipped = {'vecs':watvecs, 'scores':scores}
        fpath = '%s/%s.pdb'%(dr,batch[trg_idx])
        write_out_pdb(fpath,watvecs,scores)

def run_place_water_part(vec_starts, np_inputs, out, vcut=0.6,grid=0.5,padding=4.0):
    kernel = build_gkernel(grid=grid, r = 1.52) 
    np_o    = np.expand_dims((np_inputs[:,0,:,:,:] > 0.5),axis=1)  
    np_c    = np.expand_dims((np_inputs[:,1,:,:,:] > 0.5),axis=1) 
    np_n    = np.expand_dims((np_inputs[:,2,:,:,:] > 0.5),axis=1) 
    np_s    = np.expand_dims((np_inputs[:,3,:,:,:] > 0.5),axis=1) 
    np_prot = np.logical_or(np.logical_or(np.logical_or(np_o,np_c),np_n),np_s)
    np_null = np.logical_not(np_prot)
   
    for trg_idx in range((out.shape[0])): 
        watvecs,scores =  place_water(out[trg_idx][0], np_null[trg_idx][0], 
                                    vec_starts[trg_idx], kernel, vcut=vcut, grid =grid,
                                    padding=padding,ignore_padding=True)
        
        is_bound = [False for i in range(len(watvecs))]
        for i,vec in enumerate(watvecs):
            diff    = vec - vec_starts[trg_idx]
            n_grids = out[trg_idx][0].shape
            for j in range(3):
                if diff[j] < (padding+2.25):
                    is_bound[i] = True
                    break
                elif diff[j] > (grid*n_grids[j] - padding-2.25):
                    is_bound[i] = True
                    break
        zipped = {'vecs':watvecs, 'scores':scores,'is_bound':is_bound}
    return zipped
def test_epoch(env, epoch, batchsize=4):
    with torch.no_grad():
        run_epoch(env, epoch, train=False, build=False, batchsize=batchsize)

def train_epoch(env, epoch, batchsize=4):
    run_epoch(env, epoch, train=True, build=False, batchsize=batchsize)

def build_epoch(env, epoch, train=True, prefix='build'):
    with torch.no_grad():
        run_epoch(env, epoch, train=train, build=True, batchsize=1 , dr='./')
def build_full_epoch(env, epoch, train=True, prefix='build',tlog=False,tlog_path='cpu_time.txt'):
    with torch.no_grad():
        run_full_epoch(env, epoch, train=train, dr='./',tlog=tlog,tlog_path='cpu_time.txt')

#run_full_epoch_start
def vecs_2_idxs(vecs):
    xs_temp = []
    ys_temp = []
    zs_temp = []
    for vec_start in vecs:
        xs_temp.append(vec_start[0])
        ys_temp.append(vec_start[1])
        zs_temp.append(vec_start[2])
    xs = list(set(xs_temp))
    ys = list(set(ys_temp))
    zs = list(set(zs_temp))
    xs.sort()
    ys.sort()
    zs.sort()
    return xs,ys,zs
def run_full_epoch(env, epoch, train=True,build=False, dr=None,tlog = False, tlog_path='cpu_time.txt'):
    idxs      = env['idxs']
    net       = env['net']
    device    = env['device']
    optimizer = env['optimizer']
    n_grid    = env['n_grid']
    padding   = env['padding']
    dbloss    = env['dbloss']
    if 'use_paths' not in env.keys():
        use_paths = False
    else:
        use_paths = env['use_paths']
    running_loss = 0.0
    total_loss = 0.0
    if tlog:
        time_log = open(tlog_path,'a')
        time_start = 0
        time_end = 0
    else:
        time_log = None
        time_start = 0
        time_end = 0

    for i, idx in enumerate(idxs):
        time_start = time.time()
        wat_dict = {'vecs':[] , 'scores':[],'is_bound':[]}
        if use_paths:
            paths = env['paths_dict'][idx]
            #vec_starts: [vec_start] (current build_batch_full_paths uses only one target for batch)
            data = build_batch_full_paths(paths,grid=0.5,n_grid=n_grid,padding=padding,dbloss=dbloss) 
        else:
            #vec_starts: [vec_start] (current build_batch_full uses only one target for batch)
            data = build_batch_full(idx,grid=0.5,n_grid=n_grid,padding=padding,dbloss=dbloss) 
       
        vec_starts_part = []
        for datum in data:
            #vs: [vec_start] (current build_batch_full uses only one target for batch)
            vs, np_inputs_part,np_answers_part,np_weights_part  = datum
            vec_starts_part.append(vs[0])
        np_inputs_part = data[0][1]
        xs,ys,zs =  vecs_2_idxs(vec_starts_part)
        len_idxs = [len(xs), len(ys), len(zs)]
        #grid_out: 2: doubleloss / 1: singleloss
        n_out_ch = 0
        if dbloss:
            n_out_ch = 2
        else:
            n_out_ch = 1
        grid_out =  np.zeros( (np_inputs_part.shape[0],n_out_ch,48*len(xs)+16 , 48*len(ys)+16 , 48*len(zs)+16 )) 
        grid_prot =  np.zeros( (np_inputs_part.shape[0], np_inputs_part.shape[1],48*len(xs)+16 , 48*len(ys)+16 , 48*len(zs)+16 )) 
        wat_dict = {'vecs':[] , 'scores':[],'is_bound':[]}
        
        for data_idx, datum in enumerate(data):
            #vec_starts: [vec_start] (current build_batch_full uses only one target for batch)
            vec_starts, np_inputs,np_answers,np_weights  = datum 

            #inputs  = torch.tensor(torch.FloatTensor(np_inputs),  device=device ,requires_grad=False, dtype=torch.float32)
            #answers = torch.tensor(torch.FloatTensor(np_answers), device=device ,requires_grad=False, dtype=torch.float32)
            #weights = torch.tensor(torch.FloatTensor(np_weights), device=device ,requires_grad=False, dtype=torch.float32)
            inputs  = torch.tensor(np_inputs, dtype=torch.float32, device=device, requires_grad=False)
            answers = torch.tensor(np_answers,dtype=torch.float32, device=device, requires_grad=False)
            weights = torch.tensor(np_weights,dtype=torch.float32, device=device, requires_grad=False)
            inputs  = inputs.to(device=device) #gpu
            answers = answers.to(device=device) #gpu
            weights = weights.to(device=device) #gpu
    
            outputs = net.eval()(inputs) 
            lossf = nn.BCEWithLogitsLoss(weight=weights)
            loss = lossf(outputs, answers)

            out = torch.sigmoid(outputs).cpu().detach().numpy()
            
            #building grid_prot, grid_out
            iidxs = [xs.index(vec_starts[0][0]),ys.index(vec_starts[0][1]),zs.index(vec_starts[0][2])]
            
            g_i = [8+48*iidx     for iidx in iidxs] #grid_all_init(x,y,z)
            g_f = [8+48*(iidx+1) for iidx in iidxs] #grid_all_final(x,y,z)
            w_i = [8            for iidx in iidxs] #wgs_init(x,y,z)
            w_f = [56           for iidx in iidxs] #wgs_final(x,y,z)      
       
            #mod for boundary 
            for iidx_i , iidx in enumerate(iidxs):
                if iidx == 0:
                    g_i[iidx_i] -= 8 
                    w_i[iidx_i] -= 8 
                elif iidx == (len_idxs[iidx_i] -1):
                    g_f[iidx_i] += 8 
                    w_f[iidx_i] += 8 
        
            grid_out[:,:,g_i[0]:g_f[0], g_i[1]:g_f[1], g_i[2]:g_f[2]] = out[:,:,w_i[0]:w_f[0], w_i[1]:w_f[1], w_i[2]:w_f[2]]
            grid_prot[:,:,g_i[0]:g_f[0], g_i[1]:g_f[1], g_i[2]:g_f[2]] = np_inputs[:,:,w_i[0]:w_f[0], w_i[1]:w_f[1], w_i[2]:w_f[2]]

        grid_vs = np.array([[xs[0],ys[0],zs[0]]]) 
        vcuts = [34,38,42]
        for vcut in vcuts:
            vc_float = float(vcut)
            wat_dict = run_place_water_part(grid_vs, grid_prot, grid_out,vcut=vc_float,grid=0.5,padding=padding )
            #fpath = '%s/%s.pdb'%(dr,idx)
            fpath = '%s_scut_%2d.pdb'%(idx,vcut)
            write_out_pdb(fpath,wat_dict['vecs'],wat_dict['scores'])
        time_end = time.time()
        if tlog == True:
            tdiff = time_end - time_start
            time_log.write('%s %15.7f\n'%(idx,tdiff))
    if tlog == True:
        time_log.close()
    
def run_epoch(env, epoch, train=True,build=False, batchsize=4 ,dr=None):
    idxs      = env['idxs']
    net       = env['net']
    device    = env['device']
    optimizer = env['optimizer']
    log       = env['log']
    n_grid    = env['n_grid']
    build_r   = env['build_r']
    dbloss    = env['dbloss']
    running_loss = 0.0
    total_loss = 0.0

    if build:
        batchs = batch_idxs(idxs, batch=1, shuffle= False)
    elif train:
        batchs = batch_idxs(idxs, batch=batchsize, shuffle= True)
    else:
        batchs = batch_idxs(idxs, batch=batchsize, shuffle= False)
    
    for i, batch in enumerate(batchs):
        if build:
            vec_starts, np_inputs,np_answers,np_weights = build_batch_all(batch,padding=build_r,dbloss=dbloss) 
        elif train: 
            vec_starts, np_inputs,np_answers,np_weights = build_batch(batch,train=True,n_grid=n_grid,dbloss=dbloss)
        else:
            vec_starts, np_inputs,np_answers,np_weights = build_batch(batch,train=False,n_grid=n_grid,dbloss=dbloss)
            

        #inputs  = torch.tensor(torch.FloatTensor(np_inputs),  device=device ,requires_grad=False, dtype=torch.float32)
        #answers = torch.tensor(torch.FloatTensor(np_answers), device=device ,requires_grad=False, dtype=torch.float32)
        #weights = torch.tensor(torch.FloatTensor(np_weights), device=device ,requires_grad=False, dtype=torch.float32)
        inputs  = torch.tensor(np_inputs, dtype=torch.float32, device=device, requires_grad=False)
        answers = torch.tensor(np_answers,dtype=torch.float32, device=device, requires_grad=False)
        weights = torch.tensor(np_weights,dtype=torch.float32, device=device, requires_grad=False)
        inputs = inputs.to(device=device) #gpu
        answers = answers.to(device=device) #gpu
        weights = weights.to(device=device) #gpu
    
        if build or (not train): 
            outputs = net.eval()(inputs) 
            lossf = nn.BCEWithLogitsLoss(weight=weights)
            loss = lossf(outputs, answers)
        else:
            optimizer.zero_grad()
            outputs = net.train()(inputs) 
            lossf = nn.BCEWithLogitsLoss(weight=weights)
            loss = lossf(outputs, answers)
            loss.backward()
            optimizer.step()

        if build:
            out = torch.sigmoid(outputs).cpu().detach().numpy()
            run_place_water(vec_starts, np_inputs, out, batch, dr ,vcut=2.00,grid=0.5 )
    
        # print statistics
        running_loss += loss.item()
        total_loss   += loss.item()
        if train and (not build):
            if i % 5 == 4:
                print('T[%d, %5d] loss: %.6f' %
                      (epoch + 1, i + 1, running_loss/5.))
                log.write('T[%d, %5d] loss: %.6f\n' %
                      (epoch + 1, i + 1, running_loss/5.))
    
                running_loss = 0.0
        print (batch)

    if (not train) and (not build):
       print('E[%d, %5d] loss: %.6f' %
             (epoch + 1, i + 1, total_loss/float(len(batchs))))
       log.write('E[%d, %5d] loss: %.6f\n' %
             (epoch + 1, i + 1, total_loss/float(len(batchs))))

def training(nets,names ,n_grid=32, build_r=32., build_prefixs = [],dbloss=False,batchsize=4):
    for idx, net in enumerate(nets):
        name = names[idx]
        curr_dir =  os.path.dirname(__file__)#new 
        cwd = os.getcwd()  #new
        os.chdir(curr_dir) #new
        if not os.access(name,0):
            os.mkdir(name)
        device = torch.device("cuda:0" if torch.cuda.is_available() else "cpu")
        net.to(device)
        if torch.cuda.device_count() >1:
            net = nn.DataParallel(net)
        
        states = glob.glob('%s/cifar_state_*.bin'%(name))
        states.sort()
    
        if len(states) > 0:
            start_epoch = int(states[-1].split('/')[-1].split('.')[0].split('_')[-1])
            net.load_state_dict(torch.load(states[-1]))
        else:
            start_epoch = 0
        os.chdir(cwd) #new        
        optimizer = optim.Adam(net.parameters(), lr=0.0001)
        
        log = open('%s/cnn_net_v%d.log'%(name,idx),'a')
        env_train = {'idxs':trainidxs, 'net':net, 'device':device,'dbloss':dbloss,
                     'optimizer':optimizer, 'log':log,'n_grid':n_grid,'build_r':build_r}
        
        env_test = {'idxs':testidxs, 'net':net, 'device':device,'dbloss':dbloss,
                     'optimizer':optimizer, 'log':log,'n_grid':n_grid,'build_r':build_r}
        
        for epoch in range(start_epoch, 1000):  # loop over the dataset multiple times
            train_epoch(env_train, epoch,batchsize=batchsize)
            test_epoch(env_test, epoch)
                 
            if epoch%5 ==4:
                test_epoch(env_test, epoch)

            if epoch%50 == 49:
                #save state:
                torch.save(net.module.state_dict(), '%s/cifar_state_%05d.bin'%(name,(epoch+1)))
    
        log.close()
        del(log)
        print('BUILD_TRAIN')
        build_epoch(env_train, 1000, train=True, prefix=build_prefixs[idx])
        print('BUILD_TEST')
        build_epoch(env_test , 1000, train=False, prefix=build_prefixs[idx])

def predict_path(nets,names,in_path,out_name,n_grid=64, padding=4.0, build_prefixs = [],dbloss=False,tlog=False):
    for idx, net in enumerate(nets):
        name = names[idx]
        curr_dir =  os.path.dirname(__file__)#new 
        cwd = os.getcwd()  #new
        os.chdir(curr_dir) #new
        if not os.access(name,0):
            os.mkdir(name)
        device = torch.device("cuda:0" if torch.cuda.is_available() else "cpu")
        net.to(device)
        if torch.cuda.device_count() >1:
            net = nn.DataParallel(net)
        
        states = glob.glob('%s/cifar_state_*.bin'%(name))
        states.sort()
    
        if len(states) > 0:
            start_epoch = int(states[-1].split('/')[-1].split('.')[0].split('_')[-1])
            net.load_state_dict(torch.load(states[-1]))
        else:
            start_epoch = 0
        os.chdir(cwd) #new        
        optimizer = optim.Adam(net.parameters(), lr=0.0001)
        pd = {}
        pd[out_name] = {'pro':[in_path] , 
                         'wat':[]} 
        env = {'idxs':[out_name], 'net':net, 'device':device,'dbloss':dbloss,
                     'optimizer':optimizer, 'n_grid':n_grid,'padding':padding,
                     'use_paths':True, 'paths_dict':pd}
        build_full_epoch(env, 1000, train=False, prefix=build_prefixs[idx],tlog=tlog)

def predict_path_cpu(nets,names,in_path,out_name,n_grid=64, padding=4.0, build_prefixs = [],dbloss=False,tlog=False,tlog_path='cpu_time.log'):
    for idx, net in enumerate(nets):
        name = names[idx]
        curr_dir =  os.path.dirname(__file__)#new 
        cwd = os.getcwd()  #new
        os.chdir(curr_dir) #new
        if not os.access(name,0):
            os.mkdir(name)
        
        states = glob.glob('%s/cifar_state_*.bin'%(name))
        states.sort()
    
        if len(states) > 0:
            start_epoch = int(states[-1].split('/')[-1].split('.')[0].split('_')[-1])
            net.load_state_dict(torch.load(states[-1],map_location='cpu'))
        else:
            start_epoch = 0
        os.chdir(cwd) #new        
        optimizer = optim.Adam(net.parameters(), lr=0.0001)
        pd = {}
        pd[out_name] = {'pro':[in_path] , 
                         'wat':[]} 
        env = {'idxs':[out_name], 'net':net, 'device':'cpu','dbloss':dbloss,
                     'optimizer':optimizer, 'n_grid':n_grid,'padding':padding,
                     'use_paths':True, 'paths_dict':pd}
        build_full_epoch(env, 1000, train=False, prefix=build_prefixs[idx],tlog=tlog,tlog_path=tlog_path)

if __name__ == '__main__':
    #idxs - For training (requires GPU)
    
    #idxs is list of idx, where idx is the file name of PDB, used like below
    #   paths = {'pro':['./pdb/pdb_protein/%s.pdb'%(idx)],
    #           'wat':['./pdb/pdb_water/%s.pdb'%(idx)]}
    # pdb_protein contains protein only PDB files.
    # pdb_water   contains water only PDB files - for training
    
    trainidxs  = []
    testidxs   = []
    
    if len(sys.argv) != 3:
        print('usage: GWCNN_cpu.py [input PDB/mmCIF] [output name]')
        
    else:    
        pro_path = sys.argv[1]
        out_name = sys.argv[2]
        tlog_path ='cpu_time.log'
        nets  = [networks.Net_v4_5_auxloss()]
        names = ['networks'] 
        predict_path_cpu(nets,names,pro_path,out_name,n_grid=64, padding=4.0, build_prefixs=names,dbloss=True,tlog = False)