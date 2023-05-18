import torch
import torch.nn as nn
from torch.utils.data import Dataset, DataLoader
import numpy as np
class RandomDataset(Dataset):
    def __init__(self,size,length):
        self.len  = length
        self.data = torch.randn(length,size)
        
    def __getitem__(self, index):
        return self.data[index]
    
    def __len__(self):
        return self.len 

class Model_pro(nn.Module):
    #channels:
    # atom types 0~7
    # sybyl types (drugscore format)
    #ds_list = ['Br'   , 'C.2', 'C.3' , 'C.ar' , 
    #           'C.cat', 'Cl' , 'F'   , 'I'    , 
    #           'Met'  , 'N.3', 'N.am', 'N.ar' , 
    #           'N.pl3', 'O.2', 'O.3' , 'O.co2', 
    #           'P.3'  , 'S.3']
    def __init__(self, grid = 0.5, n_grid =[48,48,48] ):
        super(Model_pro, self).__init__()
        self.grid   = grid
        self.n_grid = n_grid

    def bbox(self, vec, b_size = 3, grid = 0.5, n_grid=[48,48,48]):
        min_pt_f = (vec-b_size)/grid
        max_pt_f = 2+(vec+b_size)/grid
        #min_pt:max_pt -> count min_pt ~ max_pt-1 -> add 1 to max_pt_f
        #int(max_pt) will be lesser than max_pt   -> add another 1  
    
        min_pt = [max(0, min(n_grid[i],int(min_pt_f[i]))) for i in range(3)]
        max_pt = [max(0, min(n_grid[i],int(max_pt_f[i]))) for i in range(3)]
        return min_pt,max_pt

    def forward(self, vecs, fgs ): #vecs contain only one ch
        pro_keys = ['C','N','O','S','P','X','M','?'] 
        #halogen(X) :size of Cl: 1.75-> 1.80 
        #Metal(M): all metal ion is same/smaller than 1.37 -> 1.40
        #unknowns(?): uses carbon size
        vdw_dic = {'C':1.7,'N':1.55,'O':1.52,'S':1.80,'P':1.95, 'X':1.80, 'M':1.40, '?':1.7} 
        output = torch.cuda.FloatTensor(8,self.n_grid[0],self.n_grid[1],self.n_grid[2]).fill_(0)
        
        for ch_idx, ch in enumerate(pro_keys):
            vdw= vdw_dic[ch]
            vec_list = vecs['pro'][ch]
            fg_list  =  fgs['pro'][ch]

            for vec_idx, vec in enumerate(vec_list):
                fg = fg_list[vec_idx]
                min_pt,max_pt =  self.bbox(vec, b_size = 1.5*vdw, grid = self.grid, n_grid=self.n_grid)
                grid = self.grid
                r = vdw
                e2 = np.exp(2)
                r2 = r*r
                #
                x = torch.arange(min_pt[0],max_pt[0],1, out=torch.cuda.LongTensor() )
                y = torch.arange(min_pt[1],max_pt[1],1, out=torch.cuda.LongTensor() )
                z = torch.arange(min_pt[2],max_pt[2],1, out=torch.cuda.LongTensor() )
                xind,yind,zind =torch.meshgrid([x,y,z])
                #
                #
                xv = grid*xind.float()
                yv = grid*yind.float()
                zv = grid*zind.float()
                #
                vec_x = torch.full_like(xv,vec[0])
                vec_y = torch.full_like(yv,vec[1])
                vec_z = torch.full_like(zv,vec[2])
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
                #
                f1  = torch.exp(-2.0*d2/r2) #short dist
                f2  = ((4.0*d2)/(e2*r2) - (12.0*d)/(e2*r) + 9./e2) #medium dist
                f3  = torch.full_like(d,0.) #long dist
                #torch.where: (contidion, true,false)
                f4  = torch.where(d<(1.5*r),f2,f3)
                mask  = torch.where(d<r,f1,f4)
                output[ch_idx, min_pt[0]:max_pt[0] , min_pt[1]:max_pt[1], min_pt[2]:max_pt[2] ] += mask[:,:,:]
                #sybyl types (-1: unknowns)
                #if fg != -1:
                #    output[fg+8, min_pt[0]:max_pt[0] , min_pt[1]:max_pt[1], min_pt[2]:max_pt[2] ] += mask[:,:,:]
        return output

class Model_wat(nn.Module):

    def __init__(self, grid = 0.5, n_grid = [48,48,48]):
        super(Model_wat, self).__init__()
        self.grid   = grid
        self.n_grid = n_grid

    def bbox(self, vec, b_size = 3, grid = 0.5, n_grid=[48,48,48]):
        min_pt_f = (vec-b_size)/grid
        max_pt_f = 2+(vec+b_size)/grid 
        #min_pt:max_pt -> count min_pt ~ max_pt-1 -> add 1 to max_pt_f
        #int(max_pt) will be lesser than max_pt   -> add another 1  
    
        min_pt = np.array( ([max(0,int(min_pt_f[i])) for i in range(3)]) )
        max_pt = np.array( ([min(n_grid[i],int(max_pt_f[i])) for i in range(3)]) )
        return min_pt,max_pt

    def forward(self, vecs,weight=25.,multi=1.5): #vecs contain only one ch
        #multi: multiplier for radius (multi=1.5 -> d<1.5*radii: positive / 
        #                              multi=1.0 -> d<1.0*radii: positice)
        output = torch.cuda.FloatTensor(2,self.n_grid[0],self.n_grid[1],self.n_grid[2]).fill_(0)
        #ouput[0]: out_wat
        #ouput[1]: grid_w
        #output = np.zeros((1,self.n_grid,self.n_grid,self.n_grid))
        
        wat_keys = ['O']
        vdw_dic = {'C':1.7,'N':1.55,'O':1.52,'S':1.80,'P':1.95, 'X':1.80, 'M':1.40, '?':1.7} 
        
        for ch_idx, ch in enumerate(wat_keys):
            vdw= vdw_dic[ch]
            vec_list = vecs['wat'][ch]

            for vec in vec_list:
                min_pt,max_pt =  self.bbox(vec, b_size = 1.5*vdw, grid = self.grid, n_grid=self.n_grid)
                grid = self.grid
                r = vdw
                e2 = np.exp(2)
                r2 = r*r
                #
                x = torch.arange(min_pt[0],max_pt[0],1, out=torch.cuda.LongTensor() )
                y = torch.arange(min_pt[1],max_pt[1],1, out=torch.cuda.LongTensor() )
                z = torch.arange(min_pt[2],max_pt[2],1, out=torch.cuda.LongTensor() )
                xind,yind,zind =torch.meshgrid([x,y,z])
                #
                #
                xv = grid*xind.float()
                yv = grid*yind.float()
                zv = grid*zind.float()
                #
                vec_x = torch.full_like(xv,vec[0])
                vec_y = torch.full_like(yv,vec[1])
                vec_z = torch.full_like(zv,vec[2])
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
                one  = torch.full_like(d,1.) #1
                zero = torch.full_like(d,0.) #0
                #
                #torch.where: (contidion, true,false)
                mask  = torch.where(d<(multi*r),one,zero)
                output[0, min_pt[0]:max_pt[0] , min_pt[1]:max_pt[1], min_pt[2]:max_pt[2] ] += mask[:,:,:]
        
        output[1,:,:,:] = output[0,:,:,:]
        on      = torch.full_like(output,1.) #1
        off     = torch.full_like(output,0.) #0
        on  [1,:,:,:] = weight
        off [1,:,:,:] = 1.    

        mask_out = torch.where(output >= 0.9,on,off)
        return mask_out
def build_cuda(vecs,fgs,grid=0.5,n_grid=48,weight=25.0,multi=1.5):
    
    #prepare model
    device = torch.device("cuda:0" if torch.cuda.is_available() else "cpu")
    model_pro = Model_pro(grid,[n_grid,n_grid,n_grid])
    model_pro.to(device)
    out_pro   =model_pro(vecs,fgs).cpu().detach().numpy()
    
    device = torch.device("cuda:0" if torch.cuda.is_available() else "cpu")
    model_wat = Model_wat(grid,[n_grid,n_grid,n_grid])
    model_wat.to(device)

    out_wat_tmp   =model_wat(vecs,weight,multi=multi).cpu().detach().numpy()

    result_pro = np.expand_dims(out_pro,axis=0)
    result_wat = np.expand_dims(np.expand_dims(out_wat_tmp[0],axis=0),axis=0)
    grid_w     = np.expand_dims(np.expand_dims(out_wat_tmp[1],axis=0),axis=0)
    
    return result_pro,result_wat,grid_w

def build_cuda_rect(vecs,fgs,grid=0.5,n_grid_arr=[48,48,48],weight=25.0,multi=1.5):
    #prepare model
    device = torch.device("cuda:0" if torch.cuda.is_available() else "cpu")
    ng = torch.tensor(n_grid_arr, device=device, requires_grad=False, dtype=torch.int32)
    ng = ng.to(device=device) #gpu
    model_pro = Model_pro(grid,ng)
    model_pro.to(device)
    out_pro   =model_pro(vecs,fgs).cpu().detach().numpy()
    
    device = torch.device("cuda:0" if torch.cuda.is_available() else "cpu")
    model_wat = Model_wat(grid,n_grid_arr)
    model_wat.to(device)

    out_wat_tmp=model_wat(vecs,weight,multi=multi).cpu().detach().numpy()
    result_pro = np.expand_dims(out_pro,axis=0)
    result_wat = np.expand_dims(np.expand_dims(out_wat_tmp[0],axis=0),axis=0)
    grid_w     = np.expand_dims(np.expand_dims(out_wat_tmp[1],axis=0),axis=0)
    
    return result_pro,result_wat,grid_w

