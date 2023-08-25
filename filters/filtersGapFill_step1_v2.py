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


class processo_gapfill(object):

    options = {
            'output_asset': 'projects/mapbiomas-workspace/AMOSTRAS/col8/CAATINGA/POS-CLASS/Gap-fill/',
            'input_asset': 'projects/mapbiomas-workspace/AMOSTRAS/col8/CAATINGA/CLASS/ClassCol8V5/',
            'input_asset': 'projects/mapbiomas-workspace/AMOSTRAS/col8/CAATINGA/POS-CLASS/misto/',
            'inputAsset7': 'projects/mapbiomas-workspace/public/collection7_1/mapbiomas_collection71_integration_v1',
            'asset_bacias_buffer' : 'projects/mapbiomas-workspace/AMOSTRAS/col7/CAATINGA/bacias_hidrograficaCaatbuffer5k',
            # 'asset_bacias_buffer' : 'projects/mapbiomas-arida/ALERTAS/auxiliar/bacias_hidrografica_caatinga',
            'classMapB' : [3, 4, 5, 9,12,13,15,18,19,20,21,22,23,24,25,26,29,30,31,32,33,36,37,38,39,40,41,42,43,44,45,46,47,48,49,50,62],
            'classNew'  : [3, 4, 3, 3,12,12,21,21,21,21,21,22,22,22,22,33,29,22,33,12,33, 21,33,33,21,21,21,21,21,21,21,21,21,21, 4,12,21]
        }


    def __init__(self, nameBacia):
        self.id_bacias = nameBacia
        self.geom_bacia = ee.FeatureCollection(self.options['asset_bacias_buffer']).filter(
                                                    ee.Filter.eq('nunivotto3', nameBacia)).first().geometry()   
        print("geometria ", len(self.geom_bacia.getInfo()['coordinates']))
        self.lstbandNames = ['classification_' + str(yy) for yy in range(1985, 2023)]
        self.years = [yy for yy in range(1985, 2023)]
        # print("lista de years \n ", self.years)

        self.version = '4'
        # self.name_imgClass = 'BACIA_' + nameBacia + '_GTB_col8'
        # self.name_imgClass = 'BACIA_' + nameBacia + '_RF_col8'
        self.name_imgClass = 'BACIA_corr_mista_' + nameBacia + '_V2'
        
        
        # https://code.earthengine.google.com/4f5c6af0912ce360a5adf69e4e6989e7
        self.imgMap7 = ee.Image(self.options['inputAsset7']).clip(self.geom_bacia)#.remap(self.options['classMapB'], 
                                                                                # self.options['classNew'])
        self.imgClass = ee.Image(self.options['input_asset'] + self.name_imgClass)
        # print(listImg.getInfo())
        #########  RECLASSIFICANDO AS CLASSES 15 E 18 PARA 21 #################  or(image2)
        # self.imgClass =ee.Image().byte()
        # for item in range(len(self.lstbandNames)):
        #     if item == 0 :
        #         # imgClass = classCol8V5.where(classCol8V5.eq(0), classCol71)
        #         image_classe = ee.Image(ee.List(listImg).get(item)).unmask(0)
        #         self.imgClass = self.imgClass.addBands(image_classe.where(
        #                                 image_classe.eq(0), self.imgMap7))
        #     else:
        #         self.imgClass = self.imgClass.addBands(ee.Image(ee.List(listImg).get(item)))
        
        # self.imgClass = self.imgClass.select(self.lstbandNames)
        print("todas as bandas \n === > ", self.imgClass.bandNames().getInfo())
        # sys.exit()
        # self.imgClass = self.imgClass.mask(self.imgClass.neq(0))  
        # o segundo processo de revisão começa na versão 3      
        
        
    def dictionary_bands(self, key, value):
        imgT = ee.Algorithms.If(
                        ee.Number(value).eq(2),
                        self.imgClass.select([key]).byte(),
                        ee.Image().rename([key]).byte().updateMask(self.imgClass.select(0))
                    )
        return ee.Image(imgT)

    def applyGapFill(self):

        lstImgMap = ee.Image().toByte()
        previousImage = None        
        for cc, yyear in enumerate(self.years):
            bandActive = 'classification_' + str(yyear)
            if yyear < 2022:                
                currentImage = self.imgClass.select(bandActive).remap(self.options['classMapB'], 
                                                        self.options['classNew']).unmask(0).rename(bandActive)
                currentMap7 = self.imgMap7.select(bandActive).remap(self.options['classMapB'], 
                                                        self.options['classNew']).rename(bandActive)
                # print("banda col 8", currentImage.bandNames().getInfo())
                # print("banda col7", currentMap7.bandNames().getInfo())

                maskGap = currentImage.eq(0)
                newBandActive = currentImage.where(maskGap, currentMap7)
                previousImage = copy.deepcopy(newBandActive)
              
            else:
                currentImage = self.imgClass.select(bandActive).remap(self.options['classMapB'], 
                                                        self.options['classNew']).unmask(0).rename(bandActive)
                
                newBandActive = currentImage.where(maskGap, previousImage)
            
            lstImgMap = lstImgMap.addBands(newBandActive)
                

        imageFilledTn = ee.Image.cat(lstImgMap).select(self.lstbandNames)
        # print("banda col 8", imageFilledTn.bandNames().getInfo())
        return imageFilledTn

    def processing_gapfill(self):


        # apply the gap fill
        imageFilled = self.applyGapFill()
        print("passou")
        # print(imageFilled.bandNames().getInfo())
        # imageFilledT0Tn = imageFilledT0Tn.set(
        #                 'version', self.version, 'biome', 'CAATINGA',
        #                 'collection', '8.0', 'id_bacia', self.id_bacias,
        #                 'sensor', 'Landsat', 'source','geodatin',
        #                 'system:footprint', self.imgClass.get('system:footprint')
        #             )
        name_toexport = 'filterGF_BACIA_'+ str(self.id_bacias) + "_V" + str(self.version)
        imageFilled = ee.Image(imageFilled).set(
                            'version', str(self.version), 
                            'biome', 'CAATINGA',
                            'source', 'geodatin',
                            'type_filter', 'gap_fill',
                            'collection', '8.0',
                            'id_bacia', self.id_bacias,
                            'sensor', 'Landsat',
                            'system:footprint' , self.imgClass.get('system:footprint')
                        )
        
        self.processoExportar(imageFilled, name_toexport)

    #exporta a imagem classificada para o asset
    def processoExportar(self, mapaRF,  nomeDesc):
        
        idasset =  self.options['output_asset'] + nomeDesc
        optExp = {
            'image': mapaRF, 
            'description': nomeDesc, 
            'assetId':idasset, 
            'region':self.geom_bacia.getInfo()['coordinates'],
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



param = {    
    'bioma': "CAATINGA", #nome do bioma setado nos metadados  
    'numeroTask': 6,
    'numeroLimit': 42,
    'conta' : {
        '0': 'caatinga01',
        '6': 'caatinga02',
        '12': 'caatinga03',
        '18': 'caatinga04',
        '24': 'caatinga05',        
        '30': 'solkan1201',
        '36': 'diegoGmail'   
    }
}

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


listaNameBacias = [
    '741','7421','7422','744','745','746','7492','751','752','753',
    '754','755','756','757','758','759','7621','7622','763','764',
    '765','766','767','771','772','773', '7741','7742','775','776',
    '777','778','76111','76116','7612',
    # '7613','7614',
    '7615','7616','7617','7618','7619'
]
cont = 30
for idbacia in listaNameBacias[:]:
    print("-----------------------------------------")
    print("----- PROCESSING BACIA {} -------".format(idbacia))

    cont = gerenciador(cont)
    aplicando_gapfill = processo_gapfill(idbacia)
    aplicando_gapfill.processing_gapfill()