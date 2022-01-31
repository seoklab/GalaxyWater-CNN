import numpy as np
from numpy import linalg
from bsphere import get_sph_pts,get_bsphere

def unitvector(alpha,beta,gamma):
    #angle between x,y : gamma
    #angle between y,z : alpha
    #angle between z,x : beta
    #fortunately, x remains same(in this distorted mess)
    a = alpha * np.pi / 180.0
    b = beta * np.pi / 180.0
    g = gamma * np.pi / 180.0

    result = np.zeros((3,3))
    #result[0] : vec x' , result[1]: vec y' , result[2]:vec z'
    result[0][0] = 1.0
    #from rotational matrix
    result[1][0] = np.cos(g)
    result[1][1] = np.sin(g)
    #from dot product
    result[2][0] = np.cos(b)
    result[2][1] = (np.cos(a) - np.cos(g)*np.cos(b))/np.sin(g)
    result[2][2] = np.sqrt(1.0 - result[2][0]**2 - result[2][1]**2)
    return result

def read_pts(fpath):
    f = open(fpath,'r')
    lines = f.readlines()
    result = []
    for line in lines:
        if not line.startswith('ATOM'):
            continue
        if not line[12:16] == ' CA ':
            continue
        vec = np.array([float(line[30+8*i:38+8*i] )for i in range(3)])
        result.append(vec)
    return result

def read_cryst(fpath):
    f = open(fpath,'r')
    lines = f.readlines()
    for line in lines:
        if line.startswith('CRYST1'):
            alpha = float(line[33:40])
            beta  = float(line[40:47])
            gamma = float(line[47:54])
            a     = float(line[6:15])
            b     = float(line[15:24])
            c     = float(line[24:33])
            abc = np.array([a,b,c])
            unit_v = unitvector(alpha,beta,gamma)
            result = { 'unit_v':unit_v,
                       'abc':abc}
            #print (np.matmul(np.transpose(unit_v),abc))
            return result

def read_mtrix(fpath):
    f = open(fpath,'r')
    lines = f.readlines()
    rotmats = []
    transvs = []
    smtry_ids = []
    
    for line in lines:
        if not line.startswith('MTRIX'):
            continue
        smtry_row = int(line[5]) -1

        smtry_id = line[7:11]
        rot_row   = [float(line[10+10*x:20+10*x]) for x in range(3)]
        trans_row = float(line[42:55])

        if smtry_id not in smtry_ids:
            smtry_ids.append(smtry_id)
            rotmats.append(np.zeros([3,3]))
            transvs.append(np.zeros([3]))
        idx = smtry_ids.index(smtry_id)
        for i in range(3):
            rotmats[idx][smtry_row][i] = rot_row[i]
        transvs[idx][smtry_row] = trans_row

    smtrylist = []
    for i in range(len(smtry_ids)):
        data = {'rot':rotmats[i],
                'tr':transvs[i]
               }
        smtrylist.append(data)
    return smtrylist

def read_smtry(fpath):
    f = open(fpath,'r')
    lines = f.readlines()
    rotmats = []
    transvs = []
    smtry_ids = []
    
    for line in lines:
        if not line.startswith('REMARK 290   SMTRY'):
            continue
        smtry_row = int(line[18]) -1
        smtry_id = line[19:23]
        rot_row   = [float(line[23+10*x:33+10*x]) for x in range(3)]
        trans_row = float(line[55:68])

        if smtry_id not in smtry_ids:
            smtry_ids.append(smtry_id)
            rotmats.append(np.zeros([3,3]))
            transvs.append(np.zeros([3]))
        idx = smtry_ids.index(smtry_id)
        for i in range(3):
            rotmats[idx][smtry_row][i] = rot_row[i]
        transvs[idx][smtry_row] = trans_row

    smtrylist = []
    for i in range(len(smtry_ids)):
        data = {'rot':rotmats[i],
                'tr':transvs[i]
               }
        smtrylist.append(data)
    return smtrylist
def unzip(x, arr,v_start):
    shape = arr.shape
    result = np.zeros_like(arr )
    prod  = arr.prod()
    for i in range(arr.shape[0]):
        prod //= arr[i]
        result[i] = x//prod
        x -= result[i]*prod
    result += v_start
    return result
def build_grid_pts(cryst, n_start=[-1,-1,-1],n_end=[1,1,1]):
    result = []
    abc    = cryst['abc']
    unit_v = cryst['unit_v']
    one = np.ones(3, dtype = np.int32)
    np_start = np.array(n_start,dtype=np.int32)
    np_end = np.array(n_end,dtype=np.int32) +one

    np_inside = np_end - np_start
    n = np.prod(np_inside)
   
    for i in range(n):
        grid_n = unzip(i,np_inside,np_start)
        grid_pt = (np.matmul(np.transpose(unit_v),(abc*grid_n)))
        result.append(grid_pt)
    return result

def write_pts(pts,outf='grid.pdb'):
    PDB_FORMAT = "HETATM%5d %4s %3s %1s%4d    %8.3f%8.3f%8.3f\n"  
    out = open(outf,'w')
    for i, pt in enumerate(pts):
        newline = PDB_FORMAT%((10000+i),' O  ','HOH','X',(9000),pt[0],pt[1],pt[2])
        out.write(newline)
    out.close()

def write_pdb_op(fpath,op,outf):
    f = open(fpath,'r')
    rot   = op['smtry']['rot']
    tr    = op['smtry']['tr']
    grid  = op['grid'] 
    lines = f.readlines()
    result = []
    out = open(outf,'w')
    for line in lines:
        if not (line.startswith('ATOM') or line.startswith('HETATM')):
            continue

        v     = np.array([float(line[30+8*i:38+8*i] )for i in range(3)])
        new_v = (np.matmul(rot,v) +tr+grid)
        newline = '%s%8.3f%8.3f%8.3f%s'%(line[:30],new_v[0],new_v[1],new_v[2],line[54:])
        out.write(newline)
    out.close()
def run_op(pts,ops):
    result = []
    for pt in pts:
        for op in ops:
            rot = op['smtry']['rot']
            tr  = op['smtry']['tr']
            grid  = op['grid']
            pt_op = (np.matmul(rot,pt) +  tr)+grid
            result.append(pt_op)
    return result

def rm_identity(ops):
    result = [] 
    for op in ops:
        rot = op['smtry']['rot']
        tr  = op['smtry']['tr'] + op['grid']
        
        mat_i = np.eye(3)
        mat_d = linalg.norm((rot-mat_i),'fro')
        tr_d  = linalg.norm(tr)
        if not(mat_d < 0.001 and tr_d < 0.001):
            result.append(op)
    return result
def find_identity(ops): 
    for op in ops:
        rot = op['smtry']['rot']
        tr  = op['smtry']['tr'] + op['grid']
        
        mat_i = np.eye(3)
        mat_d = linalg.norm((rot-mat_i),'fro')
        tr_d  = linalg.norm(tr)
        if mat_d < 0.001 and tr_d < 0.001:
            return op


def bsp_operation(bsp , smtrys, grids):
    result = [] #contains #{'center':c , 'radius':r}
    c = bsp['c']
    r = bsp['r']
    for smtry in smtrys:

        rot = smtry['rot']
        tr  = smtry['tr']
        sm_c = (np.matmul(rot,c) +  tr)

        for grid in grids:
            new_c = sm_c + grid
            newsph = {'c':new_c , 
                      'r':r}
            dat = {'sph':newsph,
                   'smtry':smtry,
                   'grid':grid}
            result.append(dat)
    return result

def do_ops_rot(ops,mat):
    #newsph = {'c':new_c , 
    #          'r':r}
    #op = {'sph':newsph,
    #       'smtry':smtry,
    #       'grid':grid}
    new_ops = []
    for op in ops:
        
        c     = op['sph']['c']
        grid  = op['grid']
        smtry = op['smtry']
        smtry_rot = op['smtry']['rot']
        smtry_tr  = op['smtry']['tr']
        new_c = np.matmul(mat,c) 
        new_g = np.matmul(mat,grid)
         
        new_sr = np.matmul(mat,smtry_rot) 
        new_st = np.matmul(mat,smtry_tr) 

        new_s ={'rot':new_sr , 'tr':new_st}
        new_sph = {'c':new_c,
                   'r':op['sph']['r']}
        data    = {'sph':new_sph,
                   'smtry':new_s,
                   'grid':new_g
                   }
        new_ops.append(data)
    return new_ops

def grid_to_box(min_pt, n_grid, grid):
    np_min_pt     = np.array(min_pt)
    np_n_grid_arr = n_grid*np.array([1.,1.,1.])
    np_max_pt = np_min_pt + grid*np_n_grid_arr 
    box = {'min':np_min_pt,
           'max':np_max_pt}
    return box

def rect_to_box(min_pt, n_grid_arr, grid):
    np_min_pt     = np.array(min_pt)
    np_n_grid_arr = np.array(n_grid_arr)
    np_max_pt = np_min_pt + grid*np_n_grid_arr 
    box = {'min':np_min_pt,
           'max':np_max_pt}
    return box
def get_box_sph(box):
    #box: { 'min':np.array([x,x,x]),
    #       'max':np.array([y,y,y])}
    #box_sph: r: np.sqrt((min-max)*(min-max))/2
    #         c: (min+max)/2
    d_v = box['max'] - box['min']
    dsq_v = d_v*d_v
    dsq = np.sum(dsq_v)
    r = np.sqrt(dsq) / 2.0
    c = 0.5*(box['max']+box['min'])
    box_sph  = {'r':r, 'c':c}
    return box_sph
def getdist(v1,v2):
    d_v = v2 - v1
    dsq_v = d_v*d_v
    dsq = np.sum(dsq_v)
    d = np.sqrt(dsq)
    return d
def find_touch_box( box , operators):

    # operator:
    # newsph = {'c':new_c , 
    #           'r':r}
    # dat = {'sph':newsph,
    #        'smtry':smtry,
    #        'grid':grid}
    op_touch = []

    for op in operators:
        is_in = sph_in_box( op['sph'], box)
        if is_in:
            op_touch.append(op)
    #box_sph = get_box_sph(box)
    #box_c = box_sph['c']
    #box_r = box_sph['r']
    #for op in operators:
    #    c= op['sph']['c']
    #    r= op['sph']['r']
    #    
    #    d = getdist(box_c,c)
    #    if d<(box_r+r):
    #        op_touch.append(op)
    return op_touch

def pt_in_line(val,min_val,max_val):
    result = True
    if val<min_val:
        result = False
    if val>max_val:
        result = False
    return result

def pt_in_box( pt, box):
    result_x = pt_in_line(pt[0],box['min'][0],box['max'][0])
    result_y = pt_in_line(pt[1],box['min'][1],box['max'][1])
    result_z = pt_in_line(pt[2],box['min'][2],box['max'][2])
    result = result_x and result_y and result_z
    return result
def find_touch_sph( sph , operators):

    op_touch = []

    for op in operators:
        is_in = sph_in_sph( op['sph'], sph)
        if is_in:
            op_touch.append(op)
    
    return op_touch

def sph_in_sph( sph1, sph2):
    #sphere touching box equals..
    #point(center of sph) touching round-edge box(roundness: radius of sphere)
    c1 = sph1['c']
    r1 = sph1['r']
    c2 = sph2['c']
    r2 = sph2['r']
    c1_pt = np.array(c1)
    c2_pt = np.array(c2)
    d = getdist(c1_pt,c2_pt)
    if d<(r1+r2):
        return True
    else:
        return False
def sph_in_box( sph, box):
    #sphere touching box equals..
    #point(center of sph) touching round-edge box(roundness: radius of sphere)
    result =False
    c = sph['c']
    r = sph['r']
    box_min = np.array(box['min'])
    box_max = np.array(box['max'])
    #1. touching with box-like area
    #1-a) min: box_min - [r,0,0] , max: box_max + [r,0,0]
    #1-b) min: box_min - [0,r,0] , max: box_max + [0,r,0]
    #1-c) min: box_min - [0,0,r] , max: box_max + [0,0,r]
    for i in range(3):
        r_v = np.zeros(3)
        r_v[i] += r
        nbox_min = box_min - r_v
        nbox_max = box_max + r_v
        nbox = {'min':nbox_min, 'max':nbox_max}
        result_tmp = pt_in_box(c,nbox)
        if result_tmp:
            #result = True
            return True
    #2. touching with round edge (total 12)
    #2-a) min_x <= c_x <=max_x, dist{ ([min_pt_y or max_pt_y],[min_pt_z or max_pt_z]) , c_yz} <=r
    #2-b) min_y <= c_y <=max_y, dist{ ([min_pt_x or max_pt_x],[min_pt_z or max_pt_z]) , c_xz} <=r
    #2-c) min_z <= c_z <=max_z, dist{ ([min_pt_x or max_pt_x],[min_pt_y or max_pt_y]) , c_xy} <=r
    for i in range(3):
        min_val = box_min[i]
        max_val = box_max[i]
        crit    = pt_in_line(c[i],min_val,max_val)
        if not crit:
            continue

        np_inside = np.array([2,2,2],dtype = np.int32)
        np_inside[i] = 1
        n = np.prod(np_inside)
        for j in range(n):
            grid_n = unzip(j,np_inside,np.zeros(3,dtype = np.int32))
            box_pt = box_min + (box_max - box_min)*grid_n
            c_pt   = np.array(c) #prevent reference copy
            box_pt[i] = 0
            c_pt[i]   = 0
            d = getdist(box_pt,c_pt)
            if d<r:
                #result = True
                return True

    #3. touching with round point (total 8)
    #dist{ ([min_pt_x or max_pt_x], [min_pt_y or max_pt_y],[min_pt_z or max_pt_z]) , c_} <=r
    np_inside = np.array([2,2,2],dtype = np.int32)
    n = np.prod(np_inside)
    for j in range(n):
        grid_n = unzip(j,np_inside,np.zeros(3, dtype=np.int32))
        box_pt = box_min + (box_max - box_min)*grid_n
        c_pt   = np.array(c) #prevent reference copy
        d = getdist(box_pt,c_pt)
        if d<r:
            #result = True
            return True
    #result = False
    return False

def get_op_pts_2(ops):
    result = []
    for op in ops:
        pt = op['sph']['c']
        result.append(pt)
    return result
def get_op_pts(ops):
    result = []
    for op in ops:
        pts = get_sph_pts(op['sph'],n=10)
        result.extend(pts)
    return result

def get_box_pts(box):
    min_v = box['min']
    max_v = box['max']

    pts= [ [min_v[0], min_v[1], min_v[2]],
           [min_v[0], min_v[1], max_v[2]],
           [min_v[0], max_v[1], min_v[2]],
           [min_v[0], max_v[1], max_v[2]],
           [max_v[0], min_v[1], min_v[2]],
           [max_v[0], min_v[1], max_v[2]],
           [max_v[0], max_v[1], min_v[2]],
           [max_v[0], max_v[1], max_v[2]] ]
    return pts

def build_all_ops(cryst_path, pts):
    cryst  = read_cryst(cryst_path) 
    smtrys = read_smtry(cryst_path)
    
    bsp     = get_bsphere(pts) #{'center':c , 'radius':r}
    grids   = build_grid_pts(cryst)
    all_op = bsp_operation(bsp , smtrys, grids)
    return all_op
if __name__ == "__main__":
    #cryst  = read_cryst('1ay7.pdb') 
    #smtrys = read_smtry('1ay7.pdb')
    #pts = read_pts('1ay7.pdb')
    
    cryst  = read_cryst('43ca.pdb') 
    smtrys = read_smtry('43ca.pdb')
    pts = read_pts('43ca.pdb')
    
    bsp =  get_bsphere(pts) #{'center':c , 'radius':r}
    box =  grid_to_box(np.array([-10.,-10.,-10.]), 104, 0.5)
    
    grids   = build_grid_pts(cryst)
    all_bsp = bsp_operation(bsp , smtrys, grids)
    touch_ops = find_touch_box( box , all_bsp)
    
    print('len_touch_ops', len(touch_ops))
    all_pts = run_op(pts,touch_ops)
    ops_pts = get_op_pts(touch_ops)
    box_pts = get_box_pts(box)
    ops_pts_2 = get_op_pts_2(all_bsp)
    
    write_pts(all_pts,outf='pts.pdb')
    write_pts(ops_pts,outf='ops.pdb')
    write_pts(ops_pts_2,outf='ops2.pdb')
    write_pts(box_pts,outf='box.pdb')
