import logging
import zipfile
import argparse
import re
import sys
import warnings
import pywt
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import sys
import base64
import zlib
import struct
import pandas as pd
import xml.etree.ElementTree as ET
plt.style.use('seaborn-pastel')
import os
from numpy import arange
from peaklist import *
import pywt




class masslist_data():
    def __init__(self, path):
        super().__init__()
        self.path = path
        self.filename= str(os.path.split(self.path)[-1])
        self.initial= self.filename[0:3]



    def __repr__(self):
        return "name of file is " +"  "+ str(os.path.split(self.path)[-1])


    def add_objects(self):
        workers = []
        if os.path.isfile(self.path):
            xx= str(os.path.split(self.path)[-1])
            print(xx)

            if xx.endswith(".mzXML"):
                print('yes')
                tree = ET.parse(self.path)
                root = tree.getroot()
                first= root.findall('.//{http://sashimi.sourceforge.net/schema_revision/mzXML_3.2}scan')
                second= root.findall('.//{http://sashimi.sourceforge.net/schema_revision/mzXML_3.2}peaks')
                for x, y in zip(first, second):
                    x= x.attrib
                    k= y.text
                    z= y.attrib
                    obj= Worker(x, k, z)
                    workers.append(obj)
        return workers





    def find_data(self):
        workers= self.add_objects()
        mlist= []
        scanlist= []
        rtlist= []
        size = []
        intensitylist= []
        filename=[]
        initial= []


        for x in workers:
            a, b = x.tag_text()
            scan, _, _, retentiontime= x.tag_dict()
            mlist.extend(a)
            scanlist.extend([scan]* len(a))
            filename.extend([self.filename]* len(a))
            initial.extend([self.initial] * len(a))
            size.append(len(a))
            rtlist.extend([retentiontime] * len(a))
            intensitylist.extend(b)
        return mlist, rtlist, scanlist, intensitylist, size, filename, initial


    def smoothing(data):
        for ii in range(2):
            (data, coeff_d) = pywt.dwt(data, 'sym3')
        return data




    def make_dataframe(self):
        a, b, c, d, e, f, g = self.find_data()
        df1= pd.DataFrame({'mz': a, 'intensity': d, 'scan': c, 'rt': b, 'filename': f, 'initial': g})
        return df1


    def peaklist_full(self):
        df1= self.make_dataframe()
        df1.to_csv(f'{self.filename}.csv')
        df1['rt']= df1['rt'].astype(float)
        minimum= df1['rt'].min()
        maximum= df1['rt'].max()
        thresold= 0.2
        ranges= arange(minimum, maximum, thresold)
        bins = pd.cut(df1['rt'], ranges)
        grouped= df1.groupby(bins, as_index= False)['intensity'].agg(np.mean)
        grouped1= df1.groupby(bins, as_index= False)['rt'].agg(np.mean)
        import plotly.express as px
        fig = px.line(x=grouped1.rt, y=grouped.intensity)
        fig.write_html(f"{self.filename}_TIC.html")
        df1['mz']= df1['mz'].astype(float)
        minimum= df1['mz'].min()
        maximum= df1['mz'].max()
        thresold= 0.005
        ranges= arange(minimum, maximum, thresold)
        bins = pd.cut(df1['mz'], ranges)
        grouped= df1.groupby(bins, as_index= False)['intensity'].agg(np.mean)
        grouped1= df1.groupby(bins, as_index= False)['mz'].agg(np.mean)
        import plotly.express as px
        fig = px.line(x=grouped1.mz, y=grouped.intensity)
        fig.write_html(f"{self.filename}_masstrace.html")
        obj_peak= Peaklist_data(df1)
        peaks= obj_peak.find_peak()
        return peaks




class Worker():
    def __init__(self, x, y, z):
        super().__init__()
        self.tag = x
        self.text = y
        self.peaktag= z

    def tag_dict(self):
        data= self.tag
        scan= int(data['num'])
        basePeakMz= data['basePeakMz']
        basePeakIntensity= data['basePeakIntensity']
        basePeakIntensity= basePeakIntensity.strip()
        basePeakIntensity= re.sub('e0', 'e+0', basePeakIntensity)
        basePeakIntensity= round(float(basePeakIntensity), 4)
        retentionTime= data['retentionTime']
        retentionTime= re.sub('[^0-9.]+', '', retentionTime)
        retentionTime= float(retentionTime)
        retentionTime= retentionTime/60
        return scan, basePeakMz, basePeakIntensity, retentionTime



    def tag_text(self):
        mz_list, intensity_list= [], []
        dd= self.peaktag
        coded= self.text
        mz_list= []
        precision = 'f'
        if dd['precision'] == 64:
            precision = 'd'

        # get endian
        endian = '!'
        if dd['byteOrder'] == 'little':
            endian = '<'
        elif dd['byteOrder'] == 'big':
            endian = '>'
        compression=None


        data = coded.encode("ascii")
        data = base64.decodebytes(data)

        mz_list= []

        if dd['compressionType'] == 'zlib':
            data = zlib.decompress(data)


        # convert from binary
        count = len(data) // struct.calcsize(endian + precision)
        data = struct.unpack(endian + precision * count, data[0:len(data)])
        points = map(list, zip(data[::2], data[1::2]))

        for x in points:

            mz_list.append(round(x[0], 4))
            intensity_list.append(round(x[-1], 4))

        mz_list= np.array(mz_list)
        intensity_list= np.array(intensity_list)
        peak= peak_loc(intensity_list)
        peak_mz= mz_list[peak]
        peak_intensity= intensity_list[peak]

        return peak_mz, peak_intensity




class Peaklist_data():

    def __init__(self, obj):
        super().__init__()
        self.mass = obj



    def find_peak(self):
        try:

            df1= self.mass
            df1= df1[df1['intensity'] > 1000]
            print('peak_picking step')
            print(df1.head())
            print(df1.info())
            df1= df1.dropna()

            df1['rt']= df1['rt'].astype(float)
            minimum= df1['rt'].min()
            print(minimum)


            maximum= df1['rt'].max()
            print(maximum)

            thresold= 0.8

            dfm= pd.DataFrame()

            ranges= arange(minimum, maximum, thresold)
            bins = pd.cut(df1['rt'], ranges)


            for c, d in df1.groupby(bins, as_index= False):
                minimum1= d['mz'].min()
                maximum1= d['mz'].max()
                thresold1= 0.005
                ranges1= arange(minimum1, maximum1, thresold1)
                bins1 = pd.cut(d['mz'], ranges1)
                grouped= d.groupby(bins1, as_index= False)

                for y, w in grouped:
                    dfx= final(w)
                    dfm= dfm.append(dfx, ignore_index= True)


            return dfm


        except Exception as e:
            pass
