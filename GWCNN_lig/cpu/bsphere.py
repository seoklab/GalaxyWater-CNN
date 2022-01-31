import numpy as np
#bounding sphere using ritter's method
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

def minmax_id(pts):
    np_pts = np.array(pts)
    min_idxs = np.argmin(np_pts,axis=0)  #[id_minx,id_miny,id_minz]
    max_idxs = np.argmax(np_pts,axis=0)  #[id_maxx,id_maxy,id_maxz]
    idxs = np.concatenate((min_idxs,max_idxs),axis=None) #[idmin,...,idmax,..]
    return idxs


def get_init_sph(pts):
    idxs = minmax_id(pts)
    maxdsq  = 0.0
    maxpair = [0,0]
    for i in range(idxs.shape[0]):
        for j in range(i):
            d_v = pts[idxs[i]] - pts[idxs[j]]
            dsq_v = d_v*d_v
            dsq = np.sum(dsq_v)
            if dsq > maxdsq:
                maxdsq = dsq
                maxpair[0] = i
                maxpair[1] = j

    radius = np.sqrt(maxdsq)/2.
    center = (pts[idxs[maxpair[0]]] +pts[idxs[maxpair[1]]])/2.
    return center,radius
    
def get_bsphere(pts):
    center,radius = get_init_sph(pts)
    c   = center
    r   = radius
    rsq = r*r
     
    for pt in pts:
        d_v = pt - c
        dsq_v = d_v*d_v
        dsq = np.sum(dsq_v)
        if dsq<rsq:
            continue
        d   = np.sqrt(dsq)
        u_v = d_v / d
        r   = (r+d) / 2.
        rsq = r*r
        c  += u_v * (d-r) /2.
    data = {'c':c, 'r':r}
    return data
def get_sph_pts(data,n=5):
    c = data['c']
    r = data['r']
    result = []
    init_v = np.array([0.,0.,r])
    for i in range(n):
        for j in range(n):
            for k in range(n):
                x1 = i /float(n)
                x2 = j /float(n)
                x3 = k /float(n)
                rotmat = get_rotmat(x1,x2,x3)
                v_tmp = np.matmul(init_v,rotmat)
                v = v_tmp+c
                result.append(v)
    return result


