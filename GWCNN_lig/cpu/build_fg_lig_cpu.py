#build grid with functional group data
# change with build_fg: build_data series gives vec_start as well 
from btorch_lig_cpu import build_cuda,build_cuda_rect
import numpy as np
import pickle
import random
from smtry import build_all_ops, find_touch_box,grid_to_box,rect_to_box,do_ops_rot
from globals import modres
pdbdir = './'
def get_rotmat(x1,x2,x3):
    #James Arvo, Fast Random Rotation Matrices (1992)
    #returns uniform distribution of rotation matrices
    # 0< x1,x2,x3 < 1
    R = np.array([ [ np.cos(2.*np.pi*x1) , np.sin(2.*np.pi*x1) , 0],
                   [-np.sin(2.*np.pi*x1) , np.cos(2.*np.pi*x1) , 0],
                   [                   0 ,                   0 , 1] ])

    v = np.array(  [ np.cos(2.*np.pi*x2)*np.sqrt(x3),
                     np.sin(2.*np.pi*x2)*np.sqrt(x3),
                     np.sqrt(1.-x3) ])[np.newaxis]
    H = np.eye(3) - 2.*np.matmul(v.T,v)

    M = -np.matmul(H,R)
    return M
def merge_all(vecs):
    #paths = {'pro':['pdb/%s_b1.pdb'%idx,'pdb/%s_b2.pdb'%idx],
    #        'wat':['water/%s_all.pdb'%idx]}
    tmp = []
    #for k in vecs.keys(): #
    for k in ['pro']:
        for ch in vecs[k].keys():
            if len(vecs[k][ch]) > 0:
                tmp.append(vecs[k][ch])
    result = np.vstack(tmp)
    return result

def build_gkernel(grid=0.5,r=1.52):
    #build gaussian_like kernel for convolution
    #of probability map and oxygen sized kernel 
    #convolve to reduce noise
    l = int( 1 + (1.5*r)/grid )
    n_grid = 2*l + 1
    kernel = np.zeros([n_grid,n_grid,n_grid])
    for v in range(n_grid**3):
        i = v % n_grid
        j = (v // n_grid)%n_grid
        k = (v // (n_grid*n_grid))%n_grid
        
        dx = grid*(i - l) 
        dy = grid*(j - l) 
        dz = grid*(k - l) 

        d =  (dx**2 + dy**2 + dz**2)**(0.5)

        kernel[i][j][k] = density_f(d,r)
    return kernel

def density_f(d,r):
    #d: euclidean distance between grid point and atom
    #r: vdw radius
    if (d <0):
        raise ValueError
    elif (d<r):
        return np.exp(-2.0*d*d/(r*r))
    elif (d <= 1.5*r): #r~1.5r, switch function to make density 0 at 1.5r
        e2= np.exp(2)
        a = (4.0*d*d) / (e2*r*r)
        b = (-12.0*d) / (e2*r)
        c = 9/(e2)
        return (a+b+c)
    else:
        return 0.0

def density_hs(d,r):
    #d: euclidean distance between grid point and atom
    #r: vdw radius
    if (d <0):
        raise ValueError
    elif (d<r):
        return 1.0
    else:
        return 0.0

def bbox(vec, b_size = 3, grid = 0.5, n_grid=48):
    min_pt_f = (vec-b_size)/grid
    max_pt_f = 2+(vec+b_size)/grid 
    #min_pt:max_pt -> count min_pt ~ max_pt-1 -> add 1 to max_pt_f
    #int(max_pt) will be lesser than max_pt   -> add another 1  

    min_pt = np.array( ([max(0,int(min_pt_f[i])) for i in range(3)]) )
    max_pt = np.array( ([min(n_grid,int(max_pt_f[i])) for i in range(3)]) )
    return min_pt,max_pt

def bbox_arr(vec, n_grid_arr,b_size = 3, grid = 0.5):
    min_pt_f = (vec-b_size)/grid
    max_pt_f = 2+(vec+b_size)/grid 
    #min_pt:max_pt -> count min_pt ~ max_pt-1 -> add 1 to max_pt_f
    #int(max_pt) will be lesser than max_pt   -> add another 1  

    min_pt = np.array( ([max(0,int(min_pt_f[i])) for i in range(3)]) )
    max_pt = np.array( ([min(n_grid_arr[i],int(max_pt_f[i])) for i in range(3)]) )
    return min_pt,max_pt
            

def slice_grid_fg(vecs,fgs,x1,x2,x3, window = 24., slice_pt= None):
    if slice_pt == None:
        slice_pt = get_slice_pt(vecs,x1,x2,x3, window)
    
    temp       = {} #
    result     = {} #
    result_fgs = {} #
    max_vdw    = 2.70 #1.8(S)*1.5

    for k in vecs.keys(): #
        result[k] = {}
        temp[k]   = {}
        result_fgs[k] = {}
        for ch in vecs[k].keys():
            result[k][ch] = []
            temp[k][ch]   = []
            result_fgs[k][ch] = []

            for idx, vec in enumerate(vecs[k][ch]):
                # change coordinate by reduce slice pt
                # for building grid
                # new vector should exist between [0,window]
                svec = vec - slice_pt
                fg   = fgs[k][ch][idx]

                if (min(svec) >= (0.0-max_vdw)) and (max(svec) <= (window+max_vdw)):
                    temp[k][ch].append(svec)
                    result_fgs[k][ch].append(fg)
    
    for k in vecs.keys(): #
        for ch in vecs[k].keys():
            if len(temp[k][ch]) == 0:
                continue
            result[k][ch] = np.vstack(temp[k][ch])
    return slice_pt, result,result_fgs
def slice_grid_full(vecs,fgs,grid=0.5,padding=16,idxx=0,slice_pt=None,n_grid=None):
    #idxx: for printing id for min,max coord only.
    #slice_pt = get_slice_pt(vecs,x1,x2,x3, window)
    slice_pt_unset = False
    if type(slice_pt) == type(None):
        if slice_pt == None:
            slice_pt_unset = True
    if type(n_grid) == type(None):
        if n_grid == None:
            slice_pt_unset = True
    if slice_pt_unset: 
        mins,maxs = get_minmax(vecs)
        med = 0.5*(mins+maxs)
        #med - padding ~ med+padding
        slice_pt      = np.array([med[i] - padding  for i in range(3)])
        n_grid = 1+ int( (2*padding - 0.01)/ grid )

    svec_max = grid*n_grid
    
    temp       = {} #
    result     = {} #
    result_fgs = {} #
    max_vdw    = 0.00

    for k in vecs.keys(): #
        result[k] = {}
        temp[k]   = {}
        result_fgs[k] = {}
        for ch in vecs[k].keys():
            result[k][ch] = []
            temp[k][ch]   = []
            result_fgs[k][ch] = []

            for idx, vec in enumerate(vecs[k][ch]):
                # change coordinate by reduce slice pt
                # for building grid
                # new vector should exist between [0,window]
                svec = vec - slice_pt
                fg   = fgs[k][ch][idx]
                lt_svec = ( svec[0] <= (svec_max+max_vdw)) and \
                          ( svec[1] <= (svec_max+max_vdw)) and \
                          ( svec[2] <= (svec_max+max_vdw))
                if (min(svec) >= (0.0-max_vdw)) and (lt_svec):
                    temp[k][ch].append(svec)
                    result_fgs[k][ch].append(fg)
    
    for k in vecs.keys(): #
        for ch in vecs[k].keys():
            if len(temp[k][ch]) == 0:
                continue
            result[k][ch] = np.vstack(temp[k][ch])
    n_grid_arr = np.array([n_grid,n_grid,n_grid])
    return slice_pt,result,result_fgs,n_grid_arr

def slice_pt_full2(vecs,grid=0.5,n_grid=64,padding=4.0):
    mins,maxs = get_minmax(vecs)
    slice_pt_range = [ [] for i in range(3)] # [[xs],[ys],[zs]]
    slice_pts = []
    for i in range(3):
        v_init = mins[i]-2*padding
        v = v_init
        while (v < (maxs[i]+2*padding)):
            slice_pt_range[i].append(v)
            v += float(n_grid) * grid - 2*padding
       
    
    for x in slice_pt_range[0]:
        for y in slice_pt_range[1]:
            for z in slice_pt_range[2]:
                slice_pts.append( np.array([x,y,z]) )

    return slice_pts

def slice_pt_full(vecs,grid=0.5,padding=16):
    
    mins,maxs = get_minmax(vecs)
    med = 0.5*(mins+maxs)
    slice_pt      = np.array([med[i] - padding  for i in range(3)])
    n_grid = 1+ int( (2*padding - 0.01)/ grid )
    return slice_pt,n_grid

def get_slice_pt(vecs,x1,x2,x3, window):
    #mins: min(all points in vecs)
    #maxs: max(all points in vecs)
    #slice_min,slice max: minimum,maximum value of lowest point in slice box
    #slice_pt: slice_min + x1(2,3) *(slice_max-slice_min)
    #slice_box: [slice_pt~slice_pt+window]
    mins,maxs = get_minmax(vecs)
    #print('mins',mins)
    #print('maxs',maxs)
    slice_min = [min(maxs[i]-window , mins[i]) for i in range(3)]
    slice_max = [max(maxs[i]-window , mins[i]) for i in range(3)]
    sv        = [x1,x2,x3]
    slice_pt = [ (slice_min[i] +sv[i]*(slice_max[i]-slice_min[i]))\
                    for i in range(3)]
    return slice_pt

def rotate(vecs,x1,x2,x3):
    #vecs:{'pro':{'C':[Nx3] , 'N':..., 'O':...},
    #       'wat':{...} }
    result = {}
    rotmat = get_rotmat(x1,x2,x3)

    for k in vecs.keys(): #
        result[k] = {}

        for ch in vecs[k].keys():
            if len(vecs[k][ch]) == 0:
                result[k][ch] = []
            else:
                result[k][ch] = np.matmul(vecs[k][ch],rotmat)
    return result

def get_minmax(vecs):
    #vecs:{'pro':{'C':[Nx3] , 'N':..., 'O':...},
    #       'wat':{...} }
    #allvec:np.array([N'x3])
    allvec = merge_all(vecs)
    mins = np.array([ min(allvec.T[i]) for i in range(3)])
    maxs = np.array([ max(allvec.T[i]) for i in range(3)])
    return mins,maxs

def total_vecs(vecs):
    #merge all vector in one list
    #for building bounding sphere
    result = []
    for k in vecs.keys():
        for k2 in vecs[k].keys():
            for v in vecs[k][k2]:
                result.append(v)
    return result
def cryst_operate(vecs,fgs,ops):
    vecs_op_tmp = {}
    vecs_op = {}
    fgs_op  = {}
    #vecs:{'pro':{'C':[Nx3] , 'N':..., 'O':...},
    #       'wat':{...} }
    
    #print(ops)
    channels = ['C','N','O','S','P','X','M','?']
    for k in vecs.keys(): #
        vecs_op_tmp[k]  = {}
        vecs_op[k]      = {}
        fgs_op[k]       = {}
        
        for ch in channels:
            vecs_op_tmp[k][ch] = []
            vecs_op[k][ch]     = []
            fgs_op[k][ch]      = []
           
            for i,v in enumerate(vecs[k][ch]):
                for op in ops:
                    rot   = op['smtry']['rot']
                    tr    = op['smtry']['tr']
                    grid  = op['grid'] 
                    new_v = (np.matmul(rot,v) +tr+grid)
                    vecs_op_tmp[k][ch].append(new_v)
                    fgs_op[k][ch].append(fgs[k][ch][i])

    for k in vecs.keys():
        for ch in channels:
            if len(vecs_op_tmp[k][ch]) == 0:
                continue
            vecs_op[k][ch] = np.vstack(vecs_op_tmp[k][ch])
    return vecs_op,fgs_op

def read_as_vecs_lig(paths):
    pdb_temp = {} #
    result   = {} #
    fgs_temp   = {} #
    fgs        = {} #
    #vecs:{'pro':{'C':[Nx3] , 'N':..., 'O':...},
    #       'wat':{...} }
    
    #channels = ['C','N','O']
    channels = ['C','N','O','S','P','X','M','?'] #changed!
    #new: 'P', 'X', 'M', '?'
    #P : phosphorus (P.3 , P.2)
    #X : unknown -> halogens
    #M : metal
    #? : unknowns
    ####ds_list = ['Br'   , 'C.2', 'C.3' , 'C.ar' , 
    ####           'C.cat', 'Cl' , 'F'   , 'I'    , 
    ####           'Met'  , 'N.3', 'N.am', 'N.ar' , 
    ####           'N.pl3', 'O.2', 'O.3' , 'O.co2', 
    ####           'P.3'  , 'S.3']
    for k in paths.keys(): #
        pdb_temp[k] = {}
        result[k]   = {}
        fgs[k]   = {}
        
        for ch in channels:
            result[k][ch]   = []
            pdb_temp[k][ch] = []
            fgs[k][ch] = []


        for fpath in paths[k]:
            ftype = fpath.split('.')[-1].strip() #pdb / mol2

            if ftype == 'pdb':
                pdb_temp_k, fgs_k = read_pdb(fpath)
                for ch in channels:
                    pdb_temp[k][ch].extend(pdb_temp_k[ch])
                    fgs[k][ch].extend(fgs_k[ch])

            elif ftype == 'mol2':
                pdb_temp_k, fgs_k = read_mol2(fpath)
                for ch in channels:
                    pdb_temp[k][ch].extend(pdb_temp_k[ch])
                    fgs[k][ch].extend(fgs_k[ch])

            elif ftype == 'cif':
                pdb_temp_k, fgs_k = read_cif(fpath)
                for ch in channels:
                    pdb_temp[k][ch].extend(pdb_temp_k[ch])
                    fgs[k][ch].extend(fgs_k[ch])

    for k in paths.keys():
        for ch in channels:
            if len(pdb_temp[k][ch]) == 0:
                continue
            result[k][ch] = np.vstack(pdb_temp[k][ch])
    return result,fgs

def read_cif(fpath):
    #_atom_site.id: mandatory
    """
    data_XXXX
    #
    loop_
    _atom_site.group_PDB (ATOM/HETATM - needed for non-compound ver)
    _atom_site.id        (mandatory)
    _atom_site.type_symbol    (atom type)
    _atom_site.label_atom_id   (atom name - needed for non-compound ver)
    _atom_site.label_alt_id    ( alt name. '.' or 'A' will be accounted.) 
    _atom_site.label_comp_id   (compound name - needed for non-compound ver)
    _atom_site.label_asym_id    (chain name)
    _atom_site.label_entity_id   (chain id -> now non-alphabet can be used)
    _atom_site.label_seq_id      (seqid)
    _atom_site.pdbx_PDB_ins_code (maybe for antibody -> ? is enough)
    _atom_site.Cartn_x            
    _atom_site.Cartn_y 
    _atom_site.Cartn_z 
    _atom_site.occupancy 
    _atom_site.B_iso_or_equiv    
    _atom_site.pdbx_formal_charge 
    _atom_site.auth_seq_id 
    _atom_site.auth_comp_id 
    _atom_site.auth_asym_id 
    _atom_site.auth_atom_id 
    _atom_site.pdbx_PDB_model_num 
    """    
    
    #ATOM   1    N  N   . ALA A 1 4   ? 22.570 -20.626 -5.602 1.00 65.75  ? 5   ALA A N   1 
    channels = ['C','N','O','S','P','X','M','?'] #changed!
    atmtypes = {}
    sybyls   = {}
    for ch in channels:
        atmtypes[ch] = []
        sybyls[ch] = [] #example: sybyl['C'] = [index of C.2, ''C.3, ...]
    f = open(fpath,'r')
    lines = f.readlines()

    loop      = [] 
    isloop    = False
    isatmline = False
    items = ['_atom_site.group_PDB', '_atom_site.id', '_atom_site.type_symbol',
             '_atom_site.label_atom_id', '_atom_site.label_alt_id', '_atom_site.label_comp_id',
             '_atom_site.label_asym_id', '_atom_site.label_entity_id', '_atom_site.label_seq_id',
             '_atom_site.pdbx_PDB_ins_code','_atom_site.Cartn_x','_atom_site.Cartn_y',
             '_atom_site.Cartn_z','_atom_site.occupancy','_atom_site.B_iso_or_equiv',
             '_atom_site.pdbx_formal_charge','_atom_site.pdbx_formal_charge','_atom_site.auth_seq_id',
             '_atom_site.auth_comp_id','_atom_site.auth_asym_id','_atom_site.auth_atom_id',
             '_atom_site.pdbx_PDB_model_num']
    idx_dict = {}         
    
    for line in lines:
        if line.startswith('data'):
            continue
        elif line.startswith('#'):
            continue    
        elif line.startswith('loop_'):
            loop      = []
            isloop    = True
            isatmline = False        
        elif line.startswith('_'):
            loop.append(line.strip())
        else:
            if isloop:
                isloop = False
                if '_atom_site.id' in loop:
                    isatmline = True
                    for i, item in enumerate(items):
                        if item in loop:
                            idx_dict[item] = i
                        else:
                            idx_dict[item] = -1
                            
                    mandatory_items = [ '_atom_site.group_PDB', '_atom_site.id',
                                        '_atom_site.type_symbol','_atom_site.label_alt_id',
                                        '_atom_site.Cartn_x','_atom_site.Cartn_y','_atom_site.Cartn_z',
                                        '_atom_site.label_comp_id','_atom_site.label_atom_id']
                    for item in mandatory_items:
                        if item not in loop:
                            isatmline = False
                            
            if isatmline: #this should not be "elif"!!!
                lsp = line.strip().split()
                

                atmgroup = lsp[ idx_dict['_atom_site.group_PDB']]
                atmno    = int( lsp[idx_dict['_atom_site.id']] )
                atmtype  = lsp[ idx_dict['_atom_site.type_symbol']]

                altname  = lsp[ idx_dict['_atom_site.label_alt_id']] 
                
                x        = float(lsp[ idx_dict['_atom_site.Cartn_x']])
                y        = float(lsp[ idx_dict['_atom_site.Cartn_y']])
                z        = float(lsp[ idx_dict['_atom_site.Cartn_z']])
                      
                resname  = lsp[ idx_dict['_atom_site.label_comp_id']] 
                atmname  = lsp[ idx_dict['_atom_site.label_atom_id']]            
                
                #chainname  = lsp[ idx_dict['_atom_site.label_asym_id']] #just in case
                #chainno    = int(lsp[ idx_dict['_atom_site.label_entity_id']]) #just in case
                #resno      = int(lsp[ idx_dict['_atom_site.label_seq_id']]) #just in case
                #res_insertion  = lsp[ idx_dict['_atom_site.pdbx_PDB_ins_code ']] #just in case
                #occupancy= float(lsp[ idx_dict['_atom_site.occupancy']]) #just in case
                #bfac     = float(lsp[ idx_dict['_atom_site.B_iso_or_equiv']]) #just in case    

                if altname not in ['.','A']:
                    continue
                if atmtype == 'H':
                    continue
                    
                vec = np.array([x,y,z])
                sybyl_idx  = get_sybyl_pdb(resname,atmname,atmtype)
               
                if atmtype in ['C','N','O','S','P']:
                    atmtypes[atmtype].append(vec) # new
                    sybyls[atmtype].append(sybyl_idx)# new
                elif atmtype in ['CA','MG','ZN','HG','MN', \
                                 'CO','FE','NI','CU','CD','CO.OH']:
                    atmtypes['M'].append(vec) # new
                    sybyls['M'].append(sybyl_idx)# new
                elif atmtype in ['F','BR','CL','I']:
                    atmtypes['X'].append(vec) # new
                    sybyls['X'].append(sybyl_idx)# new
                else:
                    atmtypes['?'].append(vec) # new
                    sybyls['?'].append(-1)# new

    return atmtypes, sybyls # -> pdb_temp[k] , fgs[k]
def read_mol2(fpath):
    channels = ['C','N','O','S','P','X','M','?'] #changed!
    atmtypes = {}
    sybyls   = {}
    
    for ch in channels:
        atmtypes[ch] = []
        sybyls[ch] = [] #example: sybyl['C'] = [index of C.2, ''C.3, ...]
    
    f = open(fpath,'r')
    lines = f.readlines()

    tp = None
    for line in lines:
        if line.startswith('@'):
            tp = line.strip().lstrip('@<TRIPOS>')
            continue
        if tp == 'ATOM':
            lsp = line.split()
            if len(lsp) != 9:
                continue
            idx     = int(lsp[0])
            atmname = lsp[5]
            atmtype = atmname.split('.')[0].upper()
            if atmname.startswith('H'):
                continue 
            vec = np.array([float(lsp[2+i]) for i in range(3)])
            sybyl_idx  = get_sybyl_lig(atmname)
           
            if atmtype in ['C','N','O','S','P']:
                atmtypes[atmtype].append(vec) # new
                sybyls[atmtype].append(sybyl_idx)# new
            elif atmtype in ['CA','MG','ZN','HG','MN', \
                             'CO','FE','NI','CU','CD','CO.OH']:
                atmtypes['M'].append(vec) # new
                sybyls['M'].append(sybyl_idx)# new
            elif atmtype in ['F','BR','CL','I']:
                atmtypes['X'].append(vec) # new
                sybyls['X'].append(sybyl_idx)# new
            else:
                atmtypes['?'].append(vec) # new
                sybyls['?'].append(-1)# new

    return atmtypes, sybyls # -> pdb_temp[k] , fgs[k]

def get_sybyl_lig(atmtype):
    #change atmtype in mol2 into drugscore form
   
     
    ds_list = ['Br'   , 'C.2', 'C.3' , 'C.ar' , 
               'C.cat', 'Cl' , 'F'   , 'I'    , 
               'Met'  , 'N.3', 'N.am', 'N.ar' , 
               'N.pl3', 'O.2', 'O.3' , 'O.co2', 
               'P.3'  , 'S.3']
    if atmtype in ds_list:
        return ds_list.index(atmtype)
    if atmtype == 'O.w':
        return ds_list.index('O.3')
    if atmtype == 'N.2':
        return ds_list.index('N.ar')
    if atmtype == 'N.4':
        return ds_list.index('N.3')
    if atmtype == 'S.2':
        return ds_list.index('S.3')
    if atmtype in ['Ca','Mg','Zn','Hg','Mn', \
                     'Co.oh','Fe','Ni','Cu','Cd']:
        return ds_list.index('Met')
    return -1

def read_pdb(fpath):
    channels = ['C','N','O','S','P','X','M','?'] #changed!
    atmtypes = {}
    sybyls   = {}
    
    for ch in channels:
        atmtypes[ch] = []
        sybyls[ch] = [] #example: sybyl['C'] = [index of C.2, ''C.3, ...]
    
     
    f = open(fpath,'r')
    lines = f.readlines()

    for line in lines:
        if not (line.startswith('ATOM') or line.startswith('HETATM')):
            continue

    
        atmtype = line[76:78].lstrip().upper()
        if atmtype == '':
            atmtype_tmp = line[12:14].strip() 
            atmtype = ''.join( [ i for i in atmtype_tmp if not i.isdigit() ] )

        resname = line[17:20]
        atmname = line[12:16]
        vec = np.array( [ float(line[30+8*i:38+8*i]) for i in range(3) ] )
        alt = line[16]
        if alt not in [' ','A']: #added
            continue             #added
        #ignore water added - 231120
        if resname in ['HOH','WAT']:
            continue
        # fix resname
        if resname in modres.keys():
            restmp = modres[resname]
            resname = restmp
        if atmtype == 'SE': #MSE
            atmname = 'S'
            atmtype = 'S'

        sybyl_idx  = get_sybyl_pdb(resname,atmname,atmtype)
      
        if atmtype in ['C','N','O','S','P']:
            atmtypes[atmtype].append(vec) # new
            sybyls[atmtype].append(sybyl_idx)# new
        elif atmtype in ['CA','MG','ZN','HG','MN', \
                         'CO','FE','NI','CU','CD']:
            atmtypes['M'].append(vec) # new
            sybyls['M'].append(sybyl_idx)# new
        elif atmtype in ['F','BR','CL','I']:
            atmtypes['X'].append(vec) # new
            sybyls['X'].append(sybyl_idx)# new
        else:
            atmtypes['?'].append(vec) # new
            sybyls['?'].append(-1)# new
            continue
    return atmtypes, sybyls # -> pdb_temp[k] , fgs[k]

def get_sybyl_pdb(resname,atmname,atmtype):
    #deal with amino acid only
    atm_st = atmname.strip()
    #0: Br -> X
    if atmtype == 'BR':
        return 0
    #1: 'C.2'
    if atm_st == 'C':
        return 1
    if resname in ['ASP','ASN']:
        if atm_st == 'CG':
            return 1
    if resname in ['GLU','GLN']:
        if atm_st == 'CD':
            return 1

    #2: 'C.3'
    if atm_st == 'CA':
        return 2
    #other sp3 C atom: check at last

    #3: 'C.ar'
    if resname in ['HIS','HID','HIE','HIP']:
        #if atm_st in ['CG','ND1','CE1','NE2','CD2']:
        if atm_st in ['CG','CE1','CD2']:
            return 3
    if resname == 'PHE':
        if atm_st in ['CG','CD1','CD2','CE1','CE2','CZ']:
            return 3
    if resname == 'TYR':
        if atm_st in ['CG','CD1','CD2','CE1','CE2']:
            return 3
    if resname == 'TRP':
        if atm_st in ['CG','CE2','CZ2','CH2','CZ3','CE3']:
            return 3 

    #4: 'C.cat':
    if resname == 'ARG':
        if atm_st == 'CZ':
            return 4

    #5: 'Cl'
    if atmtype == 'CL':
        return 5
    #6: 'F'
    if atmtype == 'F':
        return 6
    #7: 'I'
    if atmtype == 'I':
        return 7
    #8: 'Met'
    if atmtype in ['CA','MG','ZN','HG','MN', \
                     'CO','FE','NI','CU','CD']:
        return 8

    #9: 'N.3'
    if resname == 'LYS':
        if atm_st == 'NZ':
            return 9
    
    #10: 'N.am'
    if atm_st == 'N':
        return 10
    if resname in ['ASN']:
        if atm_st == 'ND2':
            return 10
    if resname in ['GLN']:
        if atm_st == 'NE2':
            return 10

    #11: 'N.ar'
    if resname in ['HIS','HID','HIE','HIP']:
        if atm_st in ['ND1','NE2']:
            return 11
    if resname in ['TRP']:
        if atm_st in ['NE1']:
            return 11
    
    #12: 'N.pl3'
    if resname in ['ARG']:
        if atm_st in ['NH1','NH2','NE']:
            return 12

    #13: 'O.2'
    if atm_st == 'O':
        return 13
    if resname in ['ASN']:
        if atm_st == 'OD1':
            return 13
    if resname in ['GLN']:
        if atm_st == 'OE1':
            return 13

    #14: 'O.3'
    if resname in ['SER']:
        if atm_st == 'OG':
            return 14
    if resname in ['THR']:
        if atm_st == 'OG1':
            return 14
    if resname in ['TYR']:
        if atm_st == 'OH':
            return 14
    
    #15: 'O.co2'
    if resname in ['ASP']:
        if atm_st in ['OD1','OD2']:
            return 15
    if resname in ['GLU']:
        if atm_st in ['OE1','OE2']:
            return 15

    if atmtype == 'P':
        return 16
    if atmtype == 'S':
        return 17
    if atmtype == 'C': #
        return 2 
    return -1 

def build_data_smtry(idx,n_grid=32,dbloss = False,lig_prob=0.5):
    #lig_prob: probability for box matched in ligand only
    paths = {'pro':['%s/pdb/pdb_protein/%s.pdb'%(pdbdir,idx),
                    '%s/lig/%s.mol2'%(pdbdir,idx)],
             'wat':['%s/pdb/pdb_water/%s.pdb'%(pdbdir,idx)]} 
    rand_val = random.random()
    if rand_val < lig_prob:
        box_paths = {'pro':['%s/lig/%s.mol2'%(pdbdir,idx)],
                 'wat':[]} 
    else:
        box_paths=paths
    cryst_path = '%s/pdb/pdb_protein/%s.pdb'%(pdbdir,idx)
    
    vecs,fgs = read_as_vecs_lig(paths) #fgs: functional group ids
    box_vecs,box_fgs = read_as_vecs_lig(box_paths) #fgs: functional group ids
    all_vec = total_vecs(vecs) 
    all_ops = build_all_ops(cryst_path, all_vec)
    #new_end
    while(True):
        i = random.random()
        j = random.random()
        k = random.random()
        i2 = random.random()
        j2 = random.random()
        k2 = random.random()
                
        rotmat =  get_rotmat(i,j,k)
        ops    =  do_ops_rot(all_ops,rotmat)
        box_vecs_rot     = rotate(box_vecs,i,j,k) 
        slice_pt         = get_slice_pt(box_vecs_rot,i2,j2,k2, window=n_grid/2.0) 
        box       = grid_to_box(np.array(slice_pt), n_grid, 1.0/2.0)
        touch_ops = find_touch_box( box , ops)
        vecs_op , fgs_op = cryst_operate(vecs,fgs,touch_ops) #TODO
        vec_start,vecs_slice,fgs_slice = slice_grid_fg(vecs_op,fgs_op,i2,j2,k2,
                                                       window=n_grid/2.0,slice_pt=slice_pt)
        
        #check the grid has more than 20 atoms (O,C,N each)
        is_avail = True
        for atm in ['O','C','N']:
            if len(vecs_slice['pro'][atm]) < 4:
                is_avail = False
                break

        if not is_avail: #how about low-prob pass? TODO
            continue

        if dbloss:
            #final answer must be channel 0!
            grid_pro,grid_wat_1,grid_w_1 = build_cuda(vecs_slice,fgs_slice,n_grid=n_grid,weight=20.,multi=1.5)
            dummy,grid_wat_2,grid_w_2 = build_cuda(vecs_slice,fgs_slice,n_grid=n_grid,weight=64.,multi=1.0)
            grid_wat = np.concatenate((grid_wat_2,grid_wat_1),axis=1)
            grid_w   = np.concatenate((grid_w_2,grid_w_1),axis=1)
            return vec_start, grid_pro,grid_wat,grid_w
        else:
            grid_pro,grid_wat,grid_w = build_cuda(vecs_slice,fgs_slice,n_grid=n_grid,weight=64.,multi=1.0)
            return vec_start, grid_pro,grid_wat,grid_w

def build_data_strict_smtry(idx,i=0.5,j=0.,k=0.,i2=0.5,j2=0.5,k2=0.5,n_grid=32,dbloss=False,lig_prob=0.5):
    
    paths = {'pro':['%s/pdb/pdb_protein/%s.pdb'%(pdbdir,idx),
                    '%s/lig/%s.mol2'%(pdbdir,idx)],
             'wat':['%s/pdb/pdb_water/%s.pdb'%(pdbdir,idx)]} 
    rand_val = random.random()
    if rand_val < lig_prob:
        box_paths = {'pro':['%s/lig/%s.mol2'%(pdbdir,idx)],
                 'wat':[]} 
    else:
        box_paths=paths
    cryst_path = '%s/pdb/pdb_protein/%s.pdb'%(pdbdir,idx)
    
    vecs,fgs = read_as_vecs_lig(paths) #fgs: functional group ids
    all_vec = total_vecs(vecs) 
    all_ops = build_all_ops(cryst_path, all_vec)
    rotmat =  get_rotmat(i,j,k)
    ops    =  do_ops_rot(all_ops,rotmat)
    
    box_vecs,box_fgs = read_as_vecs_lig(box_paths) #fgs: functional group ids
    box_vecs_rot     = rotate(box_vecs,i,j,k) 
    slice_pt         = get_slice_pt(box_vecs_rot,i2,j2,k2, window=n_grid/2.0) 

    box       =  grid_to_box(np.array(slice_pt), n_grid, 1.0/2.0)
    touch_ops = find_touch_box( box , ops)
    vecs_op , fgs_op = cryst_operate(vecs,fgs,touch_ops) #TODO
    vec_start,vecs_slice,fgs_slice = slice_grid_fg(vecs_op,fgs_op,i2,j2,k2,
                                                   window=n_grid/2.0,slice_pt=slice_pt)

    if dbloss:
        #final answer must be channel 0!
        grid_pro,grid_wat_1,grid_w_1 = build_cuda(vecs_slice,fgs_slice,n_grid=n_grid,weight=20.,multi=1.5)
        dummy,grid_wat_2,grid_w_2 = build_cuda(vecs_slice,fgs_slice,n_grid=n_grid,weight=64.,multi=1.0)
        grid_wat = np.concatenate((grid_wat_2,grid_wat_1),axis=1)
        grid_w   = np.concatenate((grid_w_2,grid_w_1),axis=1)
        return vec_start, grid_pro,grid_wat,grid_w
    else:
        grid_pro,grid_wat,grid_w = build_cuda(vecs_slice,fgs_slice,n_grid=n_grid,weight=64.,multi=1.0)
        return vec_start, grid_pro,grid_wat,grid_w


def build_data_strict_all_smtry(idx,padding=32.0,i=0.5,j=0.,k=0.,i2=0.5,j2=0.5,k2=0.5,dbloss=False,lig_prob=0.5):
    paths = {'pro':['%s/pdb/pdb_protein/%s.pdb'%(pdbdir,idx),
                    '%s/lig/%s.mol2'%(pdbdir,idx)],
             'wat':['%s/pdb/pdb_water/%s.pdb'%(pdbdir,idx)]} 
    rand_val = random.random()
    if rand_val < lig_prob:
        box_paths = {'pro':['%s/lig/%s.mol2'%(pdbdir,idx)],
                 'wat':[]} 
    else:
        box_paths=paths
    cryst_path = '%s/pdb/pdb_protein/%s.pdb'%(pdbdir,idx)
    vecs,fgs = read_as_vecs_lig(paths) #fgs: functional group ids
    all_vec = total_vecs(vecs) 
    all_ops = build_all_ops(cryst_path, all_vec)

    vecs_rot   = rotate(vecs,i,j,k) 
    slice_pt, n_grid = slice_pt_full(vecs_rot,grid=0.5,padding=16,idxx=0)
    box       =  grid_to_box(np.array(slice_pt), n_grid, 1.0/2.0)
    touch_ops = find_touch_box( box , all_ops)
    vecs_op , fgs_op = cryst_operate(vecs,fgs,touch_ops) #TODO
    vec_start,vecs_slice,fgs_slice,n_grid_arr = slice_grid_full(vecs_op,fgs_op,padding=padding,idxx=idx,slice_pt=slice_pt,n_grid=n_grid)
    
    if dbloss:
        #final answer must be channel 0!
        grid_pro,grid_wat_1,grid_w_1 = build_cuda_rect(vecs_slice,fgs_slice,n_grid_arr=n_grid_arr,weight=20.,multi=1.5)
        dummy,grid_wat_2,grid_w_2 = build_cuda_rect(vecs_slice,fgs_slice,n_grid_arr=n_grid_arr,weight=64.,multi=1.0)
        grid_wat = np.concatenate((grid_wat_2,grid_wat_1),axis=1)
        grid_w   = np.concatenate((grid_w_2,grid_w_1),axis=1)
        return vec_start, grid_pro,grid_wat,grid_w
    else:
        grid_pro,grid_wat,grid_w = build_cuda_rect(vecs_slice,fgs_slice,n_grid_arr=n_grid_arr,weight=64.,multi=1.0)
        return vec_start, grid_pro,grid_wat,grid_w

def build_data_strict_paths(paths,grid=0.5,n_grid=64,padding=4.0,dbloss=False):
    #build data for full protein
    #paths = {'pro':['%s/wkgb/%s/complex2.pdb'%(pdbdir,idx)],
    #        'wat':['%s/pdb/pdb_water/%s.pdb'%(pdbdir,idx)]}
    if 'pro' not in paths.keys():
        print ('paths : {"pro":[pro_paths...],"wat":[wat_paths(=answer,for training, optional)] }')
        raise
    if type(paths['pro']) != type([]):
        print ('paths : {"pro":[pro_paths...],"wat":[wat_paths(=answer,for training, optional)] }')
        raise
    if 'wat' not in paths.keys():
        paths['wat'] = []
    vecs,fgs = read_as_vecs_lig(paths) #fgs: functional group ids
    all_vec = total_vecs(vecs) 
    
    slice_pts = slice_pt_full2(vecs,grid=grid, n_grid=n_grid, padding=padding)
    result = []
    for slice_pt in slice_pts:
        vec_start,vecs_slice,fgs_slice,n_grid_arr = slice_grid_full(vecs,fgs,padding=padding,slice_pt=slice_pt,n_grid=n_grid)
        
        if dbloss:
            #final answer must be channel 0!
            grid_pro,grid_wat_1,grid_w_1 = build_cuda_rect(vecs_slice,fgs_slice,n_grid_arr=n_grid_arr,weight=20.,multi=1.5)
            dummy,grid_wat_2,grid_w_2 = build_cuda_rect(vecs_slice,fgs_slice,n_grid_arr=n_grid_arr,weight=64.,multi=1.0)
            grid_wat = np.concatenate((grid_wat_2,grid_wat_1),axis=1)
            grid_w   = np.concatenate((grid_w_2,grid_w_1),axis=1)
            result.append( (vec_start, grid_pro,grid_wat,grid_w) )
        else:
            grid_pro,grid_wat,grid_w = build_cuda_rect(vecs_slice,fgs_slice,n_grid_arr=n_grid_arr,weight=64.,multi=1.0)
            result.append( (vec_start, grid_pro,grid_wat,grid_w) )
    return result

def build_data_strict_all(idx,padding=32.0,i=0.5,j=0.,k=0.,i2=0.5,j2=0.5,k2=0.5,dbloss=False,lig_prob=0.5):
    paths = {'pro':['%s/%s/%s_protein_noHOH.pdb'%(pdbdir,idx,idx),
                    '%s/%s/%s_ligand.mol2'%(pdbdir,idx,idx)],
             'wat':['%s/%s/%s_full_wat.pdb'%(pdbdir,idx,idx)]} #bcut: 40
    rand_val = random.random()
    if rand_val < lig_prob:
        box_paths = {'pro':['%s/%s/%s_ligand.mol2'%(pdbdir,idx,idx)],
                 'wat':[]} 
    else:
        box_paths=paths
    vecs,fgs = read_as_vecs_lig(paths) #fgs: functional group ids
    vec_start,vecs_slice,fgs_slice,n_grid_arr = slice_grid_full(vecs,fgs,padding=padding,idxx=idx)
    
    if dbloss:
        #final answer must be channel 0!
        grid_pro,grid_wat_1,grid_w_1 = build_cuda_rect(vecs_slice,fgs_slice,n_grid_arr=n_grid_arr,weight=20.,multi=1.5)
        dummy,grid_wat_2,grid_w_2 = build_cuda_rect(vecs_slice,fgs_slice,n_grid_arr=n_grid_arr,weight=64.,multi=1.0)
        grid_wat = np.concatenate((grid_wat_2,grid_wat_1),axis=1)
        grid_w   = np.concatenate((grid_w_2,grid_w_1),axis=1)
        return vec_start, grid_pro,grid_wat,grid_w
    else:
        grid_pro,grid_wat,grid_w = build_cuda_rect(vecs_slice,fgs_slice,n_grid_arr=n_grid_arr,weight=64.,multi=1.0)
        return vec_start, grid_pro,grid_wat,grid_w


