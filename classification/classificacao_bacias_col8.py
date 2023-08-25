#!/usr/bin/env python2
# -*- coding: utf-8 -*-

'''
#SCRIPT DE CLASSIFICACAO POR BACIA
#Produzido por Geodatin - Dados e Geoinformacao
#DISTRIBUIDO COM GPLv2
'''

import ee
import os 
import gee
import json
import csv
import copy
import sys
import math
import arqParametros as arqParams 
import collections
collections.Callable = collections.abc.Callable
try:
    ee.Initialize()
    print('The Earth Engine package initialized successfully!')
except ee.EEException as e:
    print('The Earth Engine package failed to initialize!')
except:
    print("Unexpected error:", sys.exc_info()[0])
    raise
# sys.setrecursionlimit(1000000000)

#============================================================
#============== FUNCTIONS FO SPECTRAL INDEX =================
#region All functions index spectral 
# Ratio Vegetation Index
def agregateBandsIndexRATIO(img):

    ratioImg = img.expression("float(b('nir') / b('red'))")\
                            .multiply(10000).rename(['ratio']).toUint16()      

    return img.addBands(ratioImg)

# Ratio Vegetation Index
def agregateBandsIndexRVI(img):

    rviImg = img.expression("float(b('red') / b('nir'))")\
                            .multiply(10000).rename(['rvi']).toUint16()       

    return img.addBands(rviImg)

def agregateBandsIndexNDVI(img):

    ndviImg = img.expression("float(b('nir') - b('red')) / (b('nir') + b('red'))")\
                            .add(1).multiply(10000).rename(['ndvi']).toUint16()       

    return img.addBands(ndviImg)

def agregateBandsIndexWater(img):

    ndwiImg = img.expression("float(b('nir') - b('swir2')) / (b('nir') + b('swir2'))")\
                            .add(1).multiply(10000).rename(['ndwi']).toUint16()       

    return img.addBands(ndwiImg)


def AutomatedWaterExtractionIndex(img):    
    awei = img.expression(
                        "float(4 * (b('green') - b('swir2')) - (0.25 * b('nir') + 2.75 * b('swir1')))"
                    ).add(1).multiply(10000).rename("awei").toUint16()          
    
    return img.addBands(awei)

def IndiceIndicadorAgua(img):    
    iiaImg = img.expression(
                        "float((b('green') - 4 *  b('nir')) / (b('green') + 4 *  b('nir')))"
                    ).add(1).multiply(10000).rename("iia").toUint16()
    
    return img.addBands(iiaImg)

def agregateBandsIndexLAI(img):
    laiImg = img.expression(
        "float(3.618 * (b('evi') - 0.118))")\
            .add(3).multiply(10000).rename(['lai']).toUint16()

    return img.addBands(laiImg)    

def agregateBandsIndexGCVI(img):    
    gcviImgA = img.expression(
        "float(b('nir')) / (b('green')) - 1")\
            .add(1).multiply(10000).rename(['gcvi']).toUint16()        
    
    return img.addBands(gcviImgA)

# Global Environment Monitoring Index GEMI 
def agregateBandsIndexGEMI(img):    
    # "( 2 * ( NIR ^2 - RED ^2) + 1.5 * NIR + 0.5 * RED ) / ( NIR + RED + 0.5 )"
    gemiImgA = img.expression(
        "float((2 * (b('nir') * b('nir') - b('red') * b('red')) + 1.5 * b('nir') + 0.5 * b('red')) / (b('nir') + b('green') + 0.5) )")\
            .add(3).multiply(10000).rename(['gemi']).toUint16()        
    
    return img.addBands(gemiImgA)

# Chlorophyll vegetation index CVI
def agregateBandsIndexCVI(img):    
    cviImgA = img.expression(
        "float(b('nir') * (b('green') / (b('blue') * b('blue'))))")\
            .multiply(10000).rename(['cvi']).toUint16()        
    
    return img.addBands(cviImgA)

# Green leaf index  GLI
def agregateBandsIndexGLI(img):    
    gliImg = img.expression(
        "float((2 * b('green') - b('red') - b('blue')) / (2 * b('green') - b('red') - b('blue')))")\
            .add(1).multiply(10000).rename(['gli']).toUint16()        
    
    return img.addBands(gliImg)

# Shape Index  IF 
def agregateBandsIndexShapeI(img):    
    shapeImgA = img.expression(
        "float((2 * b('red') - b('green') - b('blue')) / (b('green') - b('blue')))")\
            .add(1).multiply(10000).rename(['shape']).toUint16()       
    
    return img.addBands(shapeImgA)

# Aerosol Free Vegetation Index (2100 nm) 
def agregateBandsIndexAFVI(img):    
    afviImgA = img.expression(
        "float((b('nir') - 0.5 * b('swir2')) / (b('nir') + 0.5 * b('swir2')))")\
            .add(1).multiply(10000).rename(['afvi']).toUint16()        
    
    return img.addBands(afviImgA)

# Advanced Vegetation Index 
def agregateBandsIndexAVI(img):    
    aviImgA = img.expression(
        "float((b('nir')* (1.0 - b('red')) * (b('nir') - b('red'))) ** 1/3)")\
            .add(1).multiply(10000).rename(['avi']).toUint16()        
    
    return img.addBands(aviImgA)

# Bare Soil Index 
def agregateBandsIndexBSI(img):    
    bsiImg = img.expression(
        "float(((b('swir1') - b('red')) - (b('nir') + b('blue'))) / ((b('swir1') + b('red')) + (b('nir') + b('blue'))))")\
            .add(1).multiply(10000).rename(['bsi']).toUint16()        
    
    return img.addBands(bsiImg)

# BRBA	Band Ratio for Built-up Area  
def agregateBandsIndexBRBA(img):    
    brbaImg = img.expression(
        "float(b('red') / b('swir1'))")\
            .multiply(10000).rename(['brba']).toUint16()        
    
    return img.addBands(brbaImg)

# DSWI5	Disease-Water Stress Index 5
def agregateBandsIndexDSWI5(img):    
    dswi5Img = img.expression(
        "float((b('nir') + b('green')) / (b('swir1') + b('red')))")\
            .multiply(10000).rename(['dswi5']).toUint16()        
    
    return img.addBands(dswi5Img)

# LSWI	Land Surface Water Index
def agregateBandsIndexLSWI(img):    
    lswiImg = img.expression(
        "float((b('nir') - b('swir1')) / (b('nir') + b('swir1')))")\
            .add(1).multiply(10000).rename(['lswi']).toUint16()        
    
    return img.addBands(lswiImg)

# MBI	Modified Bare Soil Index
def agregateBandsIndexMBI(img):    
    mbiImg = img.expression(
        "float(((b('swir1') - b('swir2') - b('nir')) / (b('swir1') + b('swir2') + b('nir'))) + 0.5)")\
            .add(1).multiply(10000).rename(['mbi']).toUint16()       
    
    return img.addBands(mbiImg)

# UI	Urban Index	urban
def agregateBandsIndexUI(img):    
    uiImg = img.expression(
        "float((b('swir2') - b('nir')) / (b('swir2') + b('nir')))")\
            .add(1).multiply(10000).rename(['ui']).toUint16()        
    
    return img.addBands(uiImg)

# OSAVI	Optimized Soil-Adjusted Vegetation Index
def agregateBandsIndexOSAVI(img):    
    osaviImg = img.expression(
        "float(b('nir') - b('red')) / (0.16 + b('nir') + b('red'))")\
            .add(1).multiply(10000).rename(['osavi']).toUint16()        
    
    return img.addBands(osaviImg)

# Normalized Difference Red/Green Redness Index  RI
def agregateBandsIndexRI(img):        
    riImg = img.expression(
        "float(b('nir') - b('green')) / (b('nir') + b('green'))")\
            .add(1).multiply(10000).rename(['ri']).toUint16()       
    
    return img.addBands(riImg)    

# Tasselled Cap - brightness 
def agregateBandsIndexBrightness(img):    
    tasselledCapImg = img.expression(
        "float(0.3037 * b('blue') + 0.2793 * b('green') + 0.4743 * b('red')  + 0.5585 * b('nir') + 0.5082 * b('swir1') +  0.1863 * b('swir2'))")\
            .add(1).multiply(10000).rename(['brightness']).toUint16() 
    
    return img.addBands(tasselledCapImg)

# Tasselled Cap - wetness 
def agregateBandsIndexwetness(img):    
    tasselledCapImg = img.expression(
        "float(0.1509 * b('blue') + 0.1973 * b('green') + 0.3279 * b('red')  + 0.3406 * b('nir') + 0.7112 * b('swir1') +  0.4572 * b('swir2'))")\
            .add(1).multiply(10000).rename(['wetness']).toUint16() 
    
    return img.addBands(tasselledCapImg)

# Moisture Stress Index (MSI)
def agregateBandsIndexMSI(img):    
    msiImg = img.expression(
        "float( b('nir') / b('swir1'))")\
            .add(1).multiply(10000).rename(['msi']).toUint16() 
    
    return img.addBands(msiImg)

def agregateBandsIndexGVMI(img):        
    gvmiImg = img.expression(
                    "float ((b('nir')  + 0.1) - (b('swir1') + 0.02)) / ((b('nir') + 0.1) + (b('swir1') + 0.02))" 
                ).add(3).multiply(10000).rename(['gvmi']).toUint16()     

    return img.addBands(gvmiImg) 

def agregateBandsIndexsPRI(img):        
    priImg = img.expression(
                            "float((b('green') - b('blue')) / (b('green') + b('blue')))"
                        ).rename(['pri'])   
    spriImg =   priImg.expression(
                            "float((b('pri') + 1) / 2)").add(1).multiply(10000).rename(['spri'])  

    return img.addBands(spriImg)


def agregateBandsIndexCO2Flux(img):        
    ndviImg = img.expression("float(b('nir') - b('swir2')) / (b('nir') + b('swir2'))").rename(['ndvi']) 
    
    priImg = img.expression(
                            "float((b('green') - b('blue')) / (b('green') + b('blue')))"
                        ).rename(['pri'])   
    spriImg =   priImg.expression(
                            "float((b('pri') + 1) / 2)").rename(['spri'])

    co2FluxImg = ndviImg.multiply(spriImg).rename(['co2flux'])   
    
    return img.addBands(co2FluxImg)


def agregateBandsTexturasGLCM(img):        
    img = img.toInt()                
    textura2 = img.select('nir').glcmTexture(3)  
    contrastnir = textura2.select('nir_contrast').toUint16()
    #
    textura2 = img.select('red').glcmTexture(3)  
    contrastred = textura2.select('red_contrast').toUint16()
    return  img.addBands(contrastnir).addBands(contrastred)

#endregion

#region functions apply Mosaic
def CalculateIndice(imagem):

    band_feat = [
            "ratio","rvi","ndwi","awei","iia",
            "gcvi","gemi","cvi","gli","shape","afvi",
            "avi","bsi","brba","dswi5","lswi","mbi","ui",
            "osavi","ri","brightness","wetness",
            "nir_contrast","red_contrast"
        ]
    
    # imagem em Int16 com valores inteiros ate 10000        
    # imageF = self.agregateBandsgetFractions(imagem)        
    # print(imageF.bandNames().getInfo())
    imageW = imagem.divide(10000)

    imageW = agregateBandsIndexRATIO(imageW)  #
    imageW = agregateBandsIndexRVI(imageW)    #    
    imageW = agregateBandsIndexWater(imageW)  #      
    imageW = AutomatedWaterExtractionIndex(imageW)  #      
    imageW = IndiceIndicadorAgua(imageW)    #      
    imageW = agregateBandsIndexGCVI(imageW)   #   
    imageW = agregateBandsIndexGEMI(imageW)
    imageW = agregateBandsIndexCVI(imageW) 
    imageW = agregateBandsIndexGLI(imageW) 
    imageW = agregateBandsIndexShapeI(imageW) 
    imageW = agregateBandsIndexAFVI(imageW) 
    imageW = agregateBandsIndexAVI(imageW) 
    imageW = agregateBandsIndexBSI(imageW) 
    imageW = agregateBandsIndexBRBA(imageW) 
    imageW = agregateBandsIndexDSWI5(imageW) 
    imageW = agregateBandsIndexLSWI(imageW) 
    imageW = agregateBandsIndexMBI(imageW) 
    imageW = agregateBandsIndexUI(imageW) 
    imageW = agregateBandsIndexRI(imageW) 
    imageW = agregateBandsIndexOSAVI(imageW)  #     
    imageW = agregateBandsIndexwetness(imageW)   #   
    imageW = agregateBandsIndexBrightness(imageW)  #       
    imageW = agregateBandsTexturasGLCM(imageW)     #

    return imagem.addBands(imageW).select(band_feat)

def calculate_indices_x_blocos(image):

    bnd_L = ['blue','green','red','nir','swir1','swir2']        
    # band_year = [bnd + '_median' for bnd in self.option['bnd_L']]
    band_year = [
            'blue_median','green_median','red_median',
            'nir_median','swir1_median','swir2_median'
        ]
    band_drys = [bnd + '_median_dry' for bnd in bnd_L]    
    band_wets = [bnd + '_median_wet' for bnd in bnd_L]
    band_std = [bnd + '_stdDev'for bnd in bnd_L]
    band_features = [
                "ratio","rvi","ndwi","awei","iia",
                "gcvi","gemi","cvi","gli","shape","afvi",
                "avi","bsi","brba","dswi5","lswi","mbi","ui",
                "osavi","ri","brightness","wetness",
                "nir_contrast","red_contrast"] # ,"ndfia"
    # band_features.extend(self.option['bnd_L'])        
    
    image_year = image.select(band_year)
    image_year = image_year.select(band_year, bnd_L)
    # print("imagem bandas index ")    
    # print("  ", image_year.bandNames().getInfo())
    image_year = CalculateIndice(image_year)    
    # print("imagem bandas index ")    
    # print("  ", image_year.bandNames().getInfo())
    bnd_corregida = [bnd + '_median' for bnd in band_features]
    image_year = image_year.select(band_features, bnd_corregida)
    # print("imagem bandas final median \n ", image_year.bandNames().getInfo())

    image_drys = image.select(band_drys)
    image_drys = image_drys.select(band_drys, bnd_L)
    image_drys = CalculateIndice(image_drys)
    bnd_corregida = [bnd + '_median_dry' for bnd in band_features]
    image_drys = image_drys.select(band_features, bnd_corregida)
    # print("imagem bandas final dry \n", image_drys.bandNames().getInfo())

    image_wets = image.select(band_wets)
    image_wets = image_wets.select(band_wets, bnd_L)
    image_wets = CalculateIndice(image_wets)
    bnd_corregida = [bnd + '_median_wet' for bnd in band_features]
    image_wets = image_wets.select(band_features, bnd_corregida)
    # print("imagem bandas final wet \n ", image_wets.bandNames().getInfo())   

    image_year =  image_year.addBands(image_drys).addBands(image_wets) 
    return image_year

#endregion
#============================================================



param = {
    'biomes': ['CAATINGA','CERRADO','MATAATLANTICA'],
    'bioma': "CAATINGA", #nome do bioma setado nos metadados
    'asset_bacias': "projects/mapbiomas-arida/ALERTAS/auxiliar/bacias_hidrografica_caatinga",
    'asset_bacias_buffer' : 'projects/mapbiomas-workspace/AMOSTRAS/col7/CAATINGA/bacias_hidrograficaCaatbuffer5k',
    'asset_IBGE': 'users/SEEGMapBiomas/bioma_1milhao_uf2015_250mil_IBGE_geo_v4_revisao_pampa_lagoas',
    'assetOut': 'projects/mapbiomas-workspace/AMOSTRAS/col8/CAATINGA/CLASS/ClassCol8V10/',
    'assetROIs': {'id':'projects/mapbiomas-workspace/AMOSTRAS/col8/CAATINGA/ROIs/coletaROIsv3N2cluster'},
    'assetROIsExt': {'id':'projects/mapbiomas-workspace/AMOSTRAS/col8/CAATINGA/ROIs/coletaROIsv7N2manual'},    
    'classMapB': [3, 4, 5, 9,12,13,15,18,19,20,21,22,23,24,25,26,29,30,31,32,33,36,37,38,39,40,41,42,43,44,45],
    'classNew': [3, 4, 3, 3,12,12,21,21,21,21,21,22,22,22,22,33,29,22,33,12,33, 21,33,33,21,21,21,21,21,21,21],
    'asset_mosaic': 'projects/nexgenmap/MapBiomas2/LANDSAT/BRAZIL/mosaics-2',
    'version': 5,
    'anoInicial': 1985,
    'anoFinal': 2022,
    'sufix': "_01",    
    'lsBandasMap': [],
    'numeroTask': 6,
    'numeroLimit': 42,
    'conta' : {
        '0': 'caatinga01',
        '7': 'caatinga02',
        '14': 'caatinga03',
        '21': 'caatinga04',
        '28': 'caatinga05',        
        '35': 'solkan1201',
        # '36': 'rodrigo',
        # '35': 'diegoGmail',    
    },
    'pmtRF': {
        'numberOfTrees': 165, 
        'variablesPerSplit': 15,
        'minLeafPopulation': 40,
        'bagFraction': 0.8,
        'seed': 0
    },
    # https://scikit-learn.org/stable/modules/ensemble.html#gradient-boosting
    'pmtGTB': {
        'numberOfTrees': 45, 
        'shrinkage': 0.1,         
        'samplingRate': 0.8, 
        'loss': "LeastSquares",#'Huber',#'LeastAbsoluteDeviation', 
        'seed': 0
    },
    'pmtSVM' : {
        'decisionProcedure' : 'Margin', 
        'kernelType' : 'RBF', 
        'shrinking' : True, 
        'gamma' : 0.001
    } 

}
# print(param.keys())
print("vai exportar em ", param['assetOut'])
# print(param['conta'].keys())
bandNames = [
    "swir1_stdDev_1","nir_stdDev_1","green_stdDev_1","ratio_median_dry","gli_median_wet","dswi5_median_dry",
    "ri_median","osavi_median","swir2_min","shape_median","mbi_median_dry","wetness_median_dry","green_median_texture_1",
    "iia_median_wet","slopeA_1","brba_median_dry","nir_median","lswi_median_wet","red_min","rvi_median","green_min",
    "gcvi_median_dry","shape_median_dry","cvi_median_dry","blue_median_dry","mbi_median","nir_median_dry_contrast",
    "swir2_median_wet","ui_median_wet","red_median_wet","avi_median","nir_stdDev","swir1_stdDev","red_median_dry",
    "gemi_median","osavi_median_dry","blue_median_dry_1","swir2_median_dry_1","brba_median","ratio_median",
    "gli_median_dry","blue_min_1","wetness_median","green_median_wet","blue_median_wet_1","brightness_median_wet",
    "blue_min","blue_median","red_median_contrast","swir1_min_1","evi_median","blue_stdDev_1","lswi_median_dry",
    "blue_median_wet","cvi_median","red_stdDev_1","shape_median_wet","red_median_dry_1","swir2_median_wet_1",
    "dswi5_median_wet","red_median_wet_1","afvi_median","ndwi_median","avi_median_wet","gli_median","evi_median_wet",
    "nir_median_dry","gvmi_median","cvi_median_wet","swir2_min_1","iia_median","ndwi_median_dry","green_min_1",
    "ri_median_dry","osavi_median_wet","green_median_dry","ui_median_dry","red_stdDev","nir_median_wet_1",
    "swir1_median_dry_1","red_median_1","nir_median_dry_1","swir1_median_wet","blue_stdDev","bsi_median",
    "slopeA","swir1_median","swir2_median","gvmi_median_dry","red_median","gemi_median_wet","lswi_median",
    "brightness_median_dry","awei_median_wet","nir_min","afvi_median_wet","nir_median_wet","evi_median_dry",
    "swir2_median_1","ndwi_median_wet","ratio_median_wet","swir2_stdDev","gcvi_median","ui_median","rvi_median_wet",
    "green_median_wet_1","ri_median_wet","nir_min_1","rvi_median_1","swir1_median_dry","blue_median_1","green_median_1",
    "avi_median_dry","gvmi_median_wet","wetness_median_wet","swir1_median_1","dswi5_median","swir2_stdDev_1",
    "awei_median","red_min_1","mbi_median_wet","brba_median_wet","green_stdDev","green_median_texture","swir1_min",
    "awei_median_dry","swir1_median_wet_1","gemi_median_dry","nir_median_1","red_median_dry_contrast","bsi_median_1",
    "bsi_median_2","nir_median_contrast","green_median_dry_1","afvi_median_dry","gcvi_median_wet","iia_median_dry",
    "brightness_median","green_median","swir2_median_dry"
]

bandasComuns = [
    'ratio_median_dry', 'gli_median_wet', 'dswi5_median_dry', 'ri_median', 'osavi_median', 
    'swir2_min', 'shape_median', 'mbi_median_dry', 'wetness_median_dry', 'iia_median_wet', 
    'brba_median_dry', 'nir_median', 'lswi_median_wet', 'red_min', 'slope', 'rvi_median', 
    'green_min', 'gcvi_median_dry', 'shape_median_dry', 'cvi_median_dry', 'blue_median_dry', 
    'mbi_median', 'swir2_median_wet', 'ui_median_wet', 'red_median_wet', 'avi_median', 
    'nir_stdDev', 'swir1_stdDev', 'red_median_dry', 'gemi_median', 'osavi_median_dry', 
    'brba_median', 'ratio_median', 'gli_median_dry', 'wetness_median', 'green_median_wet', 
    'brightness_median_wet', 'blue_min', 'blue_median', 'lswi_median_dry', 'blue_median_wet', 
    'cvi_median', 'shape_median_wet', 'dswi5_median_wet', 'afvi_median', 'ndwi_median', 
    'avi_median_wet', 'gli_median', 'nir_median_dry', 'cvi_median_wet', 'iia_median', 
    'ndwi_median_dry', 'ri_median_dry', 'osavi_median_wet', 'green_median_dry', 
    'ui_median_dry', 'red_stdDev', 'swir1_median_wet', 'blue_stdDev', 'bsi_median', 'swir1_median', 
    'swir2_median', 'red_median', 'gemi_median_wet', 'lswi_median', 'brightness_median_dry', 
    'awei_median_wet', 'nir_min', 'afvi_median_wet', 'nir_median_wet', 'ndwi_median_wet', 
    'ratio_median_wet', 'swir2_stdDev', 'gcvi_median', 'ui_median', 'rvi_median_wet', 
    'ri_median_wet', 'swir1_median_dry', 'avi_median_dry', 'wetness_median_wet', 'dswi5_median',
    'awei_median',  'mbi_median_wet', 'brba_median_wet', 'green_stdDev', 'green_median_texture',
    'swir1_min', 'awei_median_dry', 'gemi_median_dry', 'afvi_median_dry', 'green_median',
    'gcvi_median_wet', 'iia_median_dry', 'brightness_median','swir2_median_dry', 'class'
]

#============================================================
#========================METODOS=============================
#============================================================

def gerenciador(cont):
    #0, 18, 36, 54]
    #=====================================#
    # gerenciador de contas para controlar# 
    # processos task no gee               #
    #=====================================#
    numberofChange = [kk for kk in param['conta'].keys()]    
    if str(cont) in numberofChange:

        print("conta ativa >> {} <<".format(param['conta'][str(cont)]))        
        gee.switch_user(param['conta'][str(cont)])
        gee.init()        
        gee.tasks(n= param['numeroTask'], return_list= True)        
    
    elif cont > param['numeroLimit']:
        cont = 0
    
    cont += 1    
    return cont

#exporta a imagem classificada para o asset
def processoExportar(mapaRF, regionB, nameB):
    nomeDesc = 'BACIA_'+ str(nameB)
    idasset =  param['assetOut'] + nomeDesc
    optExp = {
        'image': mapaRF, 
        'description': nomeDesc, 
        'assetId':idasset, 
        'region':regionB.getInfo(), #['coordinates']
        'scale': 30, 
        'maxPixels': 1e13,
        "pyramidingPolicy":{".default": "mode"}
    }
    task = ee.batch.Export.image.toAsset(**optExp)
    task.start() 
    print("salvando ... " + nomeDesc + "..!")
    # print(task.status())
    for keys, vals in dict(task.status()).items():
        print ( "  {} : {}".format(keys, vals))

def process_reduce_ROIsXclass(featColROIs, classVal, threhold):
    #12': 1304, '15': 1247, '18': 1280, '22': 1635, '3': 1928, '33': 1361, '4': 1378
    lst_class = [3,4,12,15,18,22,33]
    lst_class.remove(classVal)
    nFeatColROIs = featColROIs.filter(ee.Filter.inList('class', lst_class))
    featColClass = featColROIs.filter(ee.Filter.eq('class', classVal))
    featColClass = featColClass.randomColumn('random', 0, 'normal')
    featColClass = featColClass.filter(ee.Filter.lte('random', threhold))

    return nFeatColROIs.merge(featColClass)


def GetPolygonsfromFolder(nBacias, baciabuffer, yyear):
    
    getlistPtos = ee.data.getList(param['assetROIs'])
    getlistPtosExt = ee.data.getList(param['assetROIsExt'])
    ColectionPtos = ee.FeatureCollection([])
    # print("bacias vizinhas ", nBacias)
    lstROIsAg = [
        '7617','7618','7619','7492','7622','7621',
        '7741','764','765','746','767','7615'
    ]
    for idAsset in getlistPtos:         
        path_ = idAsset.get('id')
        lsFile =  path_.split("/")
        name = lsFile[-1]
        newName = name.split('_')
        # print(newName[0])
        if str(newName[0]) in nBacias and str(newName[1]) == str(yyear):
            # print("lindo ", yyear, " <> ", name)
            FeatTemp = ee.FeatureCollection(path_).filterBounds(baciabuffer)  
            if nBacias in ['7612','7613','76116']:
                FeatTemp = process_reduce_ROIsXclass(FeatTemp, 12, 0.7)
            else:
                FeatTemp = process_reduce_ROIsXclass(FeatTemp, 12, 0.3)
            # reduce Floresta 
            FeatTemp = process_reduce_ROIsXclass(FeatTemp, 3, 0.6)
            # reduce agricultura 
            FeatTemp = process_reduce_ROIsXclass(FeatTemp, 18, 0.3) 

            # print(FeatTemp.size().getInfo())
            ColectionPtos = ColectionPtos.merge(FeatTemp.select(bandasComuns))

    # if yyear in [2016, 2021]:
    #     print("yyear ", yyear, " adicionando ")
    #     for idAsset in getlistPtosExt:         
    #         path_ = idAsset.get('id')
    #         lsFile =  path_.split("/")
    #         name = lsFile[-1]
    #         newName = name.split('_')
    #         # print(newName[0])
    #         if str(newName[0]) in nBacias:
    #             FeatTemp = ee.FeatureCollection(path_).filterBounds(baciabuffer)
    #             # print(FeatTemp.size().getInfo())
    #             ColectionPtos = ColectionPtos.merge(FeatTemp.select(bandasComuns))

    ColectionPtos = ee.FeatureCollection(ColectionPtos) 
    return  ColectionPtos

def FiltrandoROIsXimportancia(nROIs, baciasAll, nbacia):

    print("aqui  ")
    limitCaat = ee.FeatureCollection('users/CartasSol/shapes/nCaatingaBff3000')
    # selecionando todas as bacias vizinhas 
    baciasB = baciasAll.filter(ee.Filter.eq('nunivotto3', nbacia))
    # limitando pelo bioma novo com buffer
    baciasB = baciasB.geometry().buffer(2000).intersection(limitCaat.geometry())
    # filtrando todo o Rois pela área construida 
    redROIs = nROIs.filterBounds(baciasB)
    mhistogram = redROIs.aggregate_histogram('class').getInfo()
    

    ROIsEnd = ee.FeatureCollection([])
    
    roisT = ee.FeatureCollection([])
    for kk, vv in mhistogram.items():
        print("class {}: == {}".format(kk, vv))
        
        roisT = redROIs.filter(ee.Filter.eq('class', int(kk)))
        roisT =roisT.randomColumn()
        
        if int(kk) == 4:
            roisT = roisT.filter(ee.Filter.gte('random',0.5))
            # print(roisT.size().getInfo())

        elif int(kk) != 21:
            roisT = roisT.filter(ee.Filter.lte('random',0.9))
            # print(roisT.size().getInfo())

        ROIsEnd = ROIsEnd.merge(roisT)
        # roisT = None
    
    return ROIsEnd

def check_dir(file_name):
    if not os.path.exists(file_name):
        arq = open(file_name, 'w+')
        arq.close()



dictPmtroArv = {
    '35': [
            '741', '746', '753', '766', '7741', '778', 
            '7616', '7617', '7618', '7619'
    ],
    '50': [
            '7422', '745', '752', '758', '7621', 
            '776', '777',  '7612', '7615'# 
    ],
    '65':  [
            '7421','744','7492','751',
            '754','755','756','757','759','7622','763','764',
            '765','767','771','772','773', '7742','775',
            '76111','76116','7614','7613'
    ]
}



ftcol_bacias = ee.FeatureCollection(param['asset_bacias'])
imagens_mosaic = ee.ImageCollection(param['asset_mosaic']).filter(
                            ee.Filter.eq('biome', param['bioma'])
                                                        ).select(arqParams.featuresreduce)
# ftcol_baciasbuffer = ee.FeatureCollection(param['asset_bacias_buffer'])
#nome das bacias que fazem parte do bioma7619
nameBacias = arqParams.listaNameBacias
print("carregando {} bacias hidrograficas ".format(len(nameBacias)))

#lista de anos
list_anos = [k for k in range(param['anoInicial'], param['anoFinal'] + 1)]
print('lista de anos entre 1985 e 2022')
param['lsBandasMap'] = ['classification_' + str(kk) for kk in list_anos]
print(param['lsBandasMap'])
list_carta = arqParams.ls_cartas

# @mosaicos: ImageCollection com os mosaicos de Mapbiomas 
bandNames = ['awei_median_dry', 'blue_stdDev', 'brightness_median', 'cvi_median_dry',]
a_file = open("features/registroBacia_Year_FeatsSel.json", "r")
dictFeatureImp = json.load(a_file)

b_file = open("features/regBacia_Year_hiperPmtrosTuningfromROIs2Y.json", 'r')
dictHiperPmtTuning = json.load(b_file)

def iterandoXBacias( _nbacia):
    # selectBacia = ftcol_bacias.filter(ee.Filter.eq('nunivotto3', _nbacia)).first()
    # https://code.earthengine.google.com/2f8ea5070d3f081a52afbcfb7a7f9d25 
    
    baciabuffer = ee.FeatureCollection(param['asset_bacias_buffer']).filter(
                            ee.Filter.eq('nunivotto3', _nbacia)).first().geometry()
    
    lsNamesBacias = arqParams.dictBaciasViz[_nbacia]
    print("lista de Bacias vizinhas", lsNamesBacias)

    imglsClasxanos = ee.Image().byte()
    mydict = None
    pmtroClass = copy.deepcopy(param['pmtGTB'])
    # print("area ", baciabuffer.area(0.1).getInfo())
    bandas_imports = []
    for cc, ano in enumerate(range(1985, 2023)):
        
        #se o ano for 2018 utilizamos os dados de 2017 para fazer a classificacao
        bandActiva = 'classification_' + str(ano)        
        print( "banda activa: " + bandActiva)
        
        bandas_lst = dictFeatureImp[_nbacia][str(ano)][:150]
        # print(lsNamesBacias)
        ROIs_toTrain = GetPolygonsfromFolder(lsNamesBacias, baciabuffer, ano)    
        bandas_ROIs = [kk for kk in ROIs_toTrain.first().propertyNames().getInfo()]  
        # print()    
        # ROIs_toTrain  = ROIs_toTrain.filter(ee.Filter.notNull(bandasComuns))
        if param['anoInicial'] == ano : #or ano == 2021            
            # pega os dados de treinamento utilizando a geometria da bacia com buffer           
            print(" Distribuição dos pontos na bacia << {} >>".format(_nbacia))
            # print(ROIs_toTrain.first().getInfo())
            # print(" ")
            print("===  {}  ===".format(ROIs_toTrain.aggregate_histogram('class').getInfo()))            
            # ===  {'12': 1304, '15': 1247, '18': 1280, '22': 1635, '3': 1928, '33': 1361, '4': 1378}  ===
        #cria o mosaico a partir do mosaico total, cortando pelo poligono da bacia    
        colmosaicMapbiomas = imagens_mosaic.filter(ee.Filter.eq('year', ano)
                                    ).filterBounds(baciabuffer).median()

        mosaicMapbiomas = calculate_indices_x_blocos(colmosaicMapbiomas)
        mosaicMapbiomas = colmosaicMapbiomas.addBands(mosaicMapbiomas)
        # print(mosaicMapbiomas.size().getInfo())
        ################################################################
        listBandsMosaic = mosaicMapbiomas.bandNames().getInfo()
        # print(listBandsMosaic)
        # print('NUMERO DE BANDAS MOSAICO ',len(listBandsMosaic) )
        # # if param['anoInicial'] == ano:
        # #     print("bandas ativas ", listBandsMosaic)
        # # for bband in lsAllprop:
        # #     if bband not in listBandsMosaic:
        # #         print("Alerta com essa banda = ", bband)
        # print('bandas importantes ', len(bandas_lst))
        bandas_filtered = [kk for kk in bandas_lst if kk in listBandsMosaic]  # 
        bandas_imports = [kk for kk in bandas_filtered if kk in bandas_ROIs]  # 
        # print("bandas cruzadas \n ",len(bandas_imports))

        # print("        ", ROIs_toTrain.first().propertyNames().getInfo())


        ###############################################################
        # print(ROIs_toTrain.size().getInfo())
        # ROIs_toTrain_filted = ROIs_toTrain.filter(ee.Filter.notNull(bandas_imports))
        # print(ROIs_toTrain_filted.size().getInfo())
        # lsAllprop = ROIs_toTrain_filted.first().propertyNames().getInfo()
        # print('PROPERTIES FEAT = ', lsAllprop)
        #cria o classificador com as especificacoes definidas acima 
        
        classifierRF = ee.Classifier.smileRandomForest(**param['pmtRF'])\
                                    .train(ROIs_toTrain, 'class', bandas_imports)

        classifiedRF = mosaicMapbiomas.classify(classifierRF, bandActiva)
        # print("parameter loading ", dictHiperPmtTuning[_nbacia])
        # # 'numberOfTrees': 50, 
        # # 'shrinkage': 0.1,    # 
        # pmtroClass['shrinkage'] = dictHiperPmtTuning[_nbacia]['2021'][0]
        # pmtroClass['numberOfTrees'] = dictHiperPmtTuning[_nbacia]['2021'][1]
        # # print("pmtros Classifier ==> ", pmtroClass)
        # # reajusta os parametros 
        # if pmtroClass['numberOfTrees'] > 35 and _nbacia in dictPmtroArv['35']:
        #     pmtroClass['numberOfTrees'] = 35
        # elif pmtroClass['numberOfTrees'] > 50 and _nbacia in dictPmtroArv['50']:
        #     pmtroClass['numberOfTrees'] = 50
        
        # print("===="*10)
        # print("pmtros Classifier Ajustado ==> ", pmtroClass)

        # ee.Classifier.smileGradientTreeBoost(numberOfTrees, shrinkage, samplingRate, maxNodes, loss, seed)
        # classifierGTB = ee.Classifier.smileGradientTreeBoost(**pmtroClass)\
        #                             .train(ROIs_toTrain, 'class', bandas_imports)

        # classifiedGTB = mosaicMapbiomas.classify(classifierGTB, bandActiva)


        # ee.Classifier.libsvm(decisionProcedure, svmType, kernelType, shrinking, degree, gamma, coef0, cost, nu, terminationEpsilon, lossEpsilon, oneClass)
        # classifierSVM = ee.Classifier.libsvm(**param['pmtSVM'])\
        #                             .train(ROIs_toTrain, 'class', bandas_imports)

        # classifiedSVM = mosaicMapbiomas.classify(classifierSVM, bandActiva)
        # print("classificando!!!! ")

        # threeClassification  = classifiedRF.addBands(classifiedGTB).addBands(classifiedSVM)
        # threeClassification = threeClassification.reduce(ee.Reducer.mode(1))
        # threeClassification = threeClassification.rename(bandActiva)

        #se for o primeiro ano cria o dicionario e seta a variavel como
        #o resultado da primeira imagem classificada
        #print("addicionando classification bands")

        if param['anoInicial'] == ano:
            print ('entrou em 1985')
            # imglsClasxanos = copy.deepcopy(classifiedGTB)
            imglsClasxanos = copy.deepcopy(classifiedRF)            
            mydict = {
                'id_bacia': _nbacia,
                'version': param['version'],
                'biome': param['bioma'],
                'collection': '8.0',
                'sensor': 'Landsat',
                'source': 'geodatin',                
            }
            imglsClasxanos = imglsClasxanos.set(mydict)
        #se nao, adiciona a imagem como uma banda a imagem que ja existia
        else:
            # print("Adicionando o mapa do ano  ", ano)
            # imglsClasxanos = imglsClasxanos.addBands(classifiedGTB)
            imglsClasxanos = imglsClasxanos.addBands(classifiedRF)
    
    # i+=1
    # print(param['lsBandasMap'])
    # seta as propriedades na imagem classificada            
    imglsClasxanos = imglsClasxanos.select(param['lsBandasMap'])    
    imglsClasxanos = imglsClasxanos.clip(baciabuffer).set("system:footprint", baciabuffer.coordinates())
    
    # nomec = _nbacia + '_' + 'GTB_col8'
    nomec = _nbacia + '_' + 'RF_col8'
    # exporta bacia
    processoExportar(imglsClasxanos, baciabuffer.coordinates(), nomec) 
    # sys.exit()


## Revisando todos as Bacias que foram feitas 
registros_proc = "registros/lsBaciasClassifyfeitasv_1.txt"
pathFolder = os.getcwd()
path_MGRS = os.path.join(pathFolder, registros_proc)
baciasFeitas = []
check_dir(path_MGRS)

arqFeitos = open(path_MGRS, 'r')
for ii in arqFeitos.readlines():    
    ii = ii[:-1]
    # print(" => " + str(ii))
    baciasFeitas.append(ii)

arqFeitos.close()
arqFeitos = open(path_MGRS, 'a+')

# 100 arvores
nameBacias = [
    '741','7421','7422','744','745','746','7492','751','752','753',
    '754','755','756','757','758','759','7621','7622','763','764',
    '765','766','767','771','772','773', '7741','7742','775','776',
    '777','778','76111','76116','7612','7614','7615','7616','7617',
    '7618','7619', '7613'
]
# nameBacias = [
#     '766', '772', '773', '7741', '7742', '776', '777', '778', 
#     '76111', '76116', '7612', '7615', '7616', '7617', '7618', 
#     '7619', '7613'
# ]

knowMapSaved = True
listBacFalta = []
cont = 0
for _nbacia in nameBacias[21:]:
    if knowMapSaved:
        try:
            nameMap = 'BACIA_' + _nbacia + '_' + 'GTB_col8'
            imgtmp = ee.Image(param['assetOut'] + nameMap)
            print("loading ", nameMap, " ", len(imgtmp.bandNames().getInfo()), "bandas ")
        except:
            listBacFalta.append(_nbacia)
    else:
        cont = gerenciador(cont) 
        print("-------------------.kmkl-------------------------------------")
        print("--------    classificando bacia " + _nbacia + "-----------------")   
        print("--------------------------------------------------------") 
        iterandoXBacias(_nbacia) 
        arqFeitos.write(_nbacia + '\n')

arqFeitos.close()


if knowMapSaved:
    print("lista de bacias que faltam \n ",listBacFalta)
    print("total ", len(listBacFalta))