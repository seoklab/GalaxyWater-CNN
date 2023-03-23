# -*- coding: utf-8 -*-
import os,sys,glob,copy,pickle
import torch
import torchvision
import torchvision.transforms as transforms
import random
from build_fg_lig_cpu import build_data_strict_all,build_gkernel,\
                     build_data_smtry,build_data_strict_smtry,build_data_strict_paths
import numpy as np
import torch.nn as nn
import torch.nn.functional as F
import networks_lig
from scipy import ndimage #new
import torch.optim as optim
def write_out_pdb(fpath,vecs,scores,scorecut = None,pro_paths=None):
    form = 'HETATM%5d  O   WAT X%4d    %8.3f%8.3f%8.3f  1.00%6.2f\n'
    length = len(vecs)
    f = open(fpath,'w')
    dump = []
    atmno_end = 0
    #write protein
    if pro_paths != None:
        for propath in pro_paths['pro']:
            f_pro = open(propath,'r')
            lines = f_pro.readlines()
            end_ter = False
            for line in lines:
                if line.startswith("ATOM") or line.startswith("HETATM"):
                    atmno_end = int(line[6:11])
                    f.write(line)
                    end_ter = False
                elif line.startswith("TER"):
                    f.write(line)
                    end_ter = True
            if not end_ter:
                f.write("TER\n")

    for i in range(length):
        write = False
        if (scorecut == None):
            write = True
        elif scores[i] > scorecut:
            write = True
        else:
            write = False
        if write:
            newline = form%((i+atmno_end+1),(i+1001),vecs[i][0],vecs[i][1],vecs[i][2],scores[i])
            dump.append( [ scores[i] , newline ] )
    #sort!
    dump_sorted = sorted(dump, key = lambda x: -1.0*x[0])
    for i in range(min(8999,len(dump_sorted))):
        f.write(dump_sorted[i][1])
    f.write("TER\nEND\n")
    f.close()
def place_water(prob, null, vec_start, kernel, vcut=1.0, grid =0.5,padding=4.0,ignore_padding=False):
    #prob, ban: 3 dimension array
    #probability: result from net
    #null: True when atom(C,N,O, water) not exists 
    #vec_start: coord of prob[0][0][0]

    #build oxygen sized gaussian kernel
    #1.52: radius of oxygen atom (in angstrom) 
    r_grid = int( 1 + (1.5*1.52)/grid ) #radius in grid debug
    d_grid = 2*r_grid+1
    result = [] 
    scores = [] 
    maxv = 999.999
   
    prob_new = prob * null
    cv  = ndimage.convolve(prob_new, kernel, mode='reflect')

    while (maxv >= vcut):
        #1. since kernel is symmetric for i operation -> ignore about convolution eqn
        #2. reflect was used for treating sliced water at border 
        
        ind = np.unravel_index(np.argmax(cv,axis=None),cv.shape) #sth like (20,30,40)
        vec = vec_start + grid*np.array(ind)
        maxv = cv[ind]

        if (maxv < vcut):
            break
       
        result.append(vec)
        scores.append(maxv)

        for v in range(d_grid**3):
            i = v % d_grid
            j = (v // d_grid)%d_grid
            k = (v // (d_grid*d_grid))%d_grid
            
            x = min(prob.shape[0]-1, max(0, i+ind[0]-r_grid)) 
            y = min(prob.shape[1]-1, max(0, j+ind[1]-r_grid))
            z = min(prob.shape[2]-1, max(0, k+ind[2]-r_grid))
            dx = grid*(i - r_grid) 
            dy = grid*(j - r_grid) 
            dz = grid*(k - r_grid) 

            d =  (dx**2 + dy**2 + dz**2)**(0.5)
            if d < (1.5*1.52):
               cv[x][y][z] = 0.0
    if ignore_padding:
        result_new = []
        scores_new = []
        for i,vec in enumerate(result):
            diff    = vec - vec_start
            n_grids = prob.shape
            ignored = False
            for j in range(len(prob.shape)):
                if diff[j] < padding:
                    ignored = True
                    break
                elif diff[j] > (grid*n_grids[j] - padding):
                    ignored = True
                    break
            if not ignored:
                result_new.append(vec)
                scores_new.append(scores[i])
        return result_new,scores_new
    else:
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
        #zipped = {'vecs':watvecs, 'scores':scores}
        
        #fpath = '%s/%s.bin'%(dr,batch[trg_idx])
        #grid_dl = open(fpath,'wb')
        #pickle.dump(zipped,grid_dl)
        
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
        if train:
            dr = 'vec_result/train/%s/'%prefix
            #if not os.access(dr,0):
            #    os.mkdir(dr)
        else:
            dr = 'vec_result/test/%s/'%prefix
            #if not os.access(dr,0):
            #    os.mkdir(dr)
        run_epoch(env, epoch, train=train, build=True, batchsize=1 , dr=dr)
def build_full_epoch(env, epoch, train=True, prefix='build'):
    with torch.no_grad():
        if train:
            dr = 'vec_result/train/%s/'%prefix
            #if not os.access(dr,0):
            #    os.mkdir(dr)
        else:
            dr = 'vec_result/test/%s/'%prefix
            #if not os.access(dr,0):
            #    os.mkdir(dr)
        run_full_epoch(env, epoch, train=train, dr=dr)

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

def run_full_epoch(env, epoch, train=True,build=False, dr=None):
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

    for i, idx in enumerate(idxs):
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

            inputs  = torch.tensor(torch.FloatTensor(np_inputs),  device=device ,requires_grad=False, dtype=torch.float32)
            answers = torch.tensor(torch.FloatTensor(np_answers), device=device ,requires_grad=False, dtype=torch.float32)
            weights = torch.tensor(torch.FloatTensor(np_weights), device=device ,requires_grad=False, dtype=torch.float32)
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
        wat_dict = run_place_water_part(grid_vs, grid_prot, grid_out,vcut=1.00,grid=0.5,padding=padding )
        
        #fpath = '%s/%s.bin'%(dr,idx)
        #grid_dl = open(fpath,'wb')
        #pickle.dump(wat_dict,grid_dl)
        #fpath = '%s/%s.pdb'%(dr,idx)
        fpath = '%s.pdb'%(idx)
        if use_paths:
            paths = env['paths_dict'][idx]
            write_out_pdb(fpath,wat_dict['vecs'],wat_dict['scores'],pro_paths = paths)
        else:
            write_out_pdb(fpath,wat_dict['vecs'],wat_dict['scores'],pro_paths = None )
    
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
    
        inputs = torch.tensor(torch.FloatTensor(np_inputs), device=device ,requires_grad=False, dtype=torch.float32)
        answers = torch.tensor(torch.FloatTensor(np_answers), device=device ,requires_grad=False, dtype=torch.float32)
        weights = torch.tensor(torch.FloatTensor(np_weights), device=device ,requires_grad=False, dtype=torch.float32)
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
            run_place_water(vec_starts, np_inputs, out, batch, dr ,vcut=0.005,grid=0.5 )
    
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


def predict_lig(nets,names ,in_path, out_name, lig_path=None, n_grid=64, padding=4.0, build_prefixs = [],dbloss=False):
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
        
        #states = glob.glob('%s/net_v%d_cifar_state_*.bin'%(name,idx))
        states = glob.glob('%s/cifar_state_*.bin'%(name))
        states.sort()
    
        if len(states) > 0:
            start_epoch = int(states[-1].split('/')[-1].split('.')[0].split('_')[-1])
            net.load_state_dict(torch.load(states[-1]))
        else:
            start_epoch = 0
        ####
        os.chdir(cwd) #new
        optimizer = optim.Adam(net.parameters(), lr=0.0001)
        pd = {}
        if lig_path == None:
            pd[out_name] = {'pro':[in_path] , 'wat':[]} 
        else:
            pd[out_name] = {'pro':[in_path,lig_path] , 'wat':[]} 
        env  = {'idxs':[out_name],  'net':net, 'device':device,'dbloss':dbloss,
                     'optimizer':optimizer, 'n_grid':n_grid,'padding':padding,
                     'use_paths':True, 'paths_dict':pd}
        build_full_epoch(env, 1000, train=False, prefix=build_prefixs[idx])

def predict_lig_cpu(nets,names ,in_path, out_name, lig_path=None, n_grid=64, padding=4.0, build_prefixs = [],dbloss=False):
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
        if lig_path == None:
            pd[out_name] = {'pro':[in_path] , 'wat':[]} 
        else:
            pd[out_name] = {'pro':[in_path,lig_path] , 'wat':[]} 
        env  = {'idxs':[out_name],  'net':net, 'device':'cpu','dbloss':dbloss,
                     'optimizer':optimizer, 'n_grid':n_grid,'padding':padding,
                     'use_paths':True, 'paths_dict':pd}
        build_full_epoch(env, 1000, train=False, prefix=build_prefixs[idx])

if __name__ == '__main__':
    #for training
    #idxs - For training (requires GPU)
    
    #idxs is list of idx, where idx is the file name of PDB, used like below
    #paths = {'pro':['./pdb/pdb_protein/%s.pdb'%(idx),
    #                './lig/%s.mol2'%(idx)],
    #         'wat':['./pdb/pdb_water/%s.pdb'%(idx)]} 
    # pdb_protein contains protein only PDB files.
    # pdb_water   contains water only PDB files - for training
    # lig         contains water only ligand mol2 files - ligands in PDB will be included too
    #                                                     necessary for training stage 
    trainidxs  = []
    testidxs   = []
    #=================================================================================
    #setting result folders
    #if not os.access('vec_result/',0):
    #    os.mkdir('vec_result/')
    #if not os.access('vec_result/train/',0):
    #    os.mkdir('vec_result/train/')
    #if not os.access('vec_result/test/',0):
    #    os.mkdir('vec_result/test/')
    #=================================================================================
    if len(sys.argv) < 3:
        print('usage: GWCNN_lig_cpu.py [input PDB/mmCIF] (input ligand mol2) [output name] ') 
        raise ValueError
        
    elif len(sys.argv) == 3:        
        pro_path = sys.argv[1]
        out_name = sys.argv[2]
        lig_path = None
        
    elif len(sys.argv) == 4:
        pro_path = sys.argv[1]
        lig_path = sys.argv[2]
        out_name = sys.argv[3]


    nets  = [networks_lig.Net_v4_5_auxloss()]
    names = ['lig_nofg']
    
    #training(nets,names,n_grid=32, build_r=16., build_prefixs=names,dbloss=True,batchsize=4)
    predict_lig_cpu(nets,names,pro_path,out_name, lig_path, n_grid=64, padding=4.0, build_prefixs=names,dbloss=True)
    
    
