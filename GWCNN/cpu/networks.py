import os,sys,glob
import torch
import torchvision
import torchvision.transforms as transforms
import random
import numpy as np
import torch.nn as nn
import torch.nn.functional as F
class Net_v4_5_auxloss(nn.Module):
    #first 10 layers: predict water (coarse)
    #second 10 layers: predict water from protein atoms & coarse water map 
    #final answer must be channel 0 !!!!! (due to place_water uses [trg,0,:,:,:])
    #initial loss only used for calculating auxillary loss 
    def __init__(self):
        super(Net_v4_5_auxloss, self).__init__()
        self.conv_i1  = nn.Conv3d(16,   64, 1,padding=0)
        self.conv1   = atr_res_12(64)
        self.conv2   = atr_res_12(64)
        self.conv3   = atr_res_12(64)
        self.conv4   = atr_res_12(64)
        self.conv5   = atr_res_12(64)
        self.conv6   = atr_res_12(64)
        self.conv7   = atr_res_12(64)
        self.conv8   = atr_res_12(64)
        self.conv9   = atr_res_12(64)
        self.conv10  = atr_res_12(64)
        self.conv_e1 = nn.Conv3d(64, 1, 1) #fc layer
        
        self.conv11   = atr_res_12(64)
        self.conv12   = atr_res_12(64)
        self.conv13   = atr_res_12(64)
        self.conv14   = atr_res_12(64)
        self.conv15   = atr_res_12(64)
        self.conv16   = atr_res_12(64)
        self.conv17   = atr_res_12(64)
        self.conv18   = atr_res_12(64)
        self.conv19   = atr_res_12(64)
        self.conv20  = atr_res_12(64)
        self.conv_e2 = nn.Conv3d(64, 1, 1) #fc layer
        
    def forward(self, x_init):
        #print('init_size',x.size())
        x = F.elu(self.conv_i1(x_init))
        x = self.conv1(x)
        x = self.conv2(x)
        x = self.conv3(x)
        x = self.conv4(x)
        x = self.conv5(x)
        x = self.conv6(x)
        x = self.conv7(x)
        x = self.conv8(x)
        x = self.conv9(x)
        x = self.conv10(x)
        coarse = self.conv_e1(x) #auxloss only

        x = self.conv11(x)
        x = self.conv12(x)
        x = self.conv13(x)
        x = self.conv14(x)
        x = self.conv15(x)
        x = self.conv16(x)
        x = self.conv17(x)
        x = self.conv18(x)
        x = self.conv19(x)
        x = self.conv20(x)
        fine = self.conv_e2(x)

        z = torch.cat((fine,coarse),1)
        return z

class atr_res_12(nn.Module):
    '''atrous_conv, dilation 1, dilation 2 w/resnet'''
    def __init__(self, ch):
        super(atr_res_12, self).__init__()
        self.conv11 = nn.Conv3d(ch,   ch, 3,padding=1,dilation=1)
        self.conv12 = nn.Conv3d(ch,   ch, 3,padding=2,dilation=2)
        self.conv13 = nn.Conv3d(2*ch, ch, 1,padding=0)
        self.bn1 = nn.BatchNorm3d(ch)

    def forward(self, x):
        x0 = x
        y1 = F.elu(self.conv11(x))
        y2 = F.elu(self.conv12(x))
        x  = F.elu(x0 + self.bn1(self.conv13(torch.cat((y1,y2),1))))
        
        return x

