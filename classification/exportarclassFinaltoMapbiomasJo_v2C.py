#!/usr/bin/env python2
# -*- coding: utf-8 -*-

##########################################################
## CRIPT DE EXPORTAÇÃO DO RESULTADO FINAL PARA O ASSET  ##
## DE mAPBIOMAS                                         ##
## Produzido por Geodatin - Dados e Geoinformação       ##
##  DISTRIBUIDO COM GPLv2                               ##
#########################################################

import ee 
import gee
import json
import csv
import sys
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
sys.setrecursionlimit(1000000000)


def gerenciador(cont, paramet):
    #0, 18, 36, 54]
    #=====================================#
    # gerenciador de contas para controlar# 
    # processos task no gee               #
    #=====================================#
    numberofChange = [kk for kk in paramet['conta'].keys()]
    
    if str(cont) in numberofChange:

        print("conta ativa >> {} <<".format(paramet['conta'][str(cont)]))        
        gee.switch_user(paramet['conta'][str(cont)])
        gee.init()        
        gee.tasks(n= paramet['numeroTask'], return_list= True)        
    
    elif cont > paramet['numeroLimit']:
        cont = 0
    
    cont += 1    
    return cont


param = {
    'input_solo': 'projects/mapbiomas-workspace/AMOSTRAS/col6/CAATINGA/classificacoes/camada_solo',
    'otherinput': 'projects/mapbiomas-workspace/AMOSTRAS/col7/CAATINGA/classification_Col71_v1',
    'inputAsset': 'projects/mapbiomas-workspace/AMOSTRAS/col8/CAATINGA/POS-CLASS/Temporal',
    'outputAsset': 'projects/mapbiomas-workspace/AMOSTRAS/col8/CAATINGA/POS-CLASS/toExp',  
    'classMapB': [3, 4, 5, 9,12,13,15,18,19,20,21,22,23,24,25,26,29,30,31,32,33,36,37,38,39,40,41,42,43,44,45],
    'classNew': [3, 4, 3, 3,12,12,21,21,21,21,21,22,22,22,22,33,29,22,33,12,33, 21,33,33,21,21,21,21,21,21,21], 
    'biome': 'CAATINGA', #configure como null se for tema transversal
    'version': "4",
    'versionInput': "12",
    'collection': 8,
    'source': 'geodatin',
    'theme': None, 
    'numeroTask': 0,
    'numeroLimit': 38,
    'conta' : {
        '0': 'caatinga01',
        '7': 'caatinga02',
        '14': 'caatinga03',
        '21': 'caatinga04',
        # '27': 'caatinga05',        
        '27': 'solkan1201',  
        # '28': 'rodrigo',
        '33': 'diegoGmail'
    }
}
mapSolo = ee.Image(param['input_solo'])
metadados = {}
# lst_bacias = ['7612','7613','76116']
bioma250mil = ee.FeatureCollection('users/CartasSol/shapes/nCaatingaBff3000').geometry()
imgColClass = ee.ImageCollection(param['inputAsset']).filter(
                        ee.Filter.eq('version', param['versionInput']))# .max()

print("Numero de imagens ", imgColClass.size().getInfo())
imgColClass = imgColClass.max()

for ii, year in enumerate(range(1985, 2023)):  #
    
    gerenciador(ii , param)
    bandaAct = 'classification_' + str(year) 
    
    # print("Banda activa: " + bandaAct)
    img_banda = 'CAATINGA-' + str(year) +  '-' + param['version']

    if year < 2021:
        camadaSolo = mapSolo.select(bandaAct)

    imgYear = imgColClass.select(bandaAct)#.remap(param['classMapB'],param['classNew'])
    imgYear = imgYear.where(camadaSolo.eq(1), camadaSolo.multiply(22))               
    imgYear = ee.Image(imgYear).clip(bioma250mil).set('biome', param['biome'])\
                    .set('year', year)\
                    .set('version', param['version'])\
                    .set('collection', param['collection'])\
                    .set('source', param['source'])\
                    .set('system:footprint', bioma250mil)    

    
    name = param['biome'] + '-' + str(year) + '-' + param['version']

    optExp = {   
        'image': imgYear.byte(), 
        'description': name, 
        'assetId': param['outputAsset'] + '/' + name, 
        'region': bioma250mil.getInfo()['coordinates'], #
        'scale': 30, 
        'maxPixels': 1e13,
        "pyramidingPolicy": {".default": "mode"}
    }

    task = ee.batch.Export.image.toAsset(**optExp)
    task.start() 
    print("salvando ... banda  " + name + "..!")
    # sys.exit()