import ee
import numpy as np
import copy




dict_trehold_OLD = {
    'std' : {
        'blueBef': 759,
        'greenBef': 846,
        'redBef': 1008,
        'NIRBef': 1733,
        'blueAft': 982,
        'greenAft': 985,
        'redAft': 1133,
        'NIRAft': 1963
    },
    'means' : {
        'blueBef': 1087,
        'greenBef': 1400,
        'redBef': 1768,
        'NIRBef': 4163,
        'blueAft': 1963,
        'greenAft': 1499,
        'redAft': 1994,
        'NIRAft': 4207
    }   
}
dict_trehold_new = {
    'std' : {
        'blueBef': 126,
        'greenBef': 174,
        'redBef': 276,
        'NIRBef': 364,
        'blueAft': 136,
        'greenAft': 187,
        'redAft': 299,
        'NIRAft': 381
    },
    'means' : {
        'blueBef': 383,
        'greenBef': 616,
        'redBef': 625,
        'NIRBef': 2435,
        'blueAft': 399,
        'greenAft': 636,
        'redAft': 661,
        'NIRAft': 2491
    },
   ,
    'quartil': {
        # qaurtil  [25%, 50%, 75%]
        'blueBef': [230, 350, 536],
        'greenBef': [426, 576, 807],
        'redBef': [288, 547, 1030],
        'NIRBef': [1970, 2414, 3084],
        'blueAft': [230, 362, 577],
        'greenAft': [430, 591, 849],
        'redAft': [282, 576, 1134],
        'NIRAft': [2004, 2478, 3124]
    }   
}



def agregarIndices (image):        
    
    bndInd = ['Blue','Green','Red','NIR','ndvi','osavi','isoil','lai','ndwi']  ##    'lai','ndwi', 'evi',

    # bandsS2 = ['Green','Red','NIR']
    
    eviImg = image.expression(
            "(2.5 * float(b('NIR') - b('Red'))/(b('NIR')  + 2.4 * b('Red') + 1))").rename('evi')
    
    laiImg = eviImg.expression(
            "(3.618 * float(b('evi') - 0.118))").rename(['lai'])  # .divide(10)
    
    ndviImg = image.normalizedDifference(['NIR', 'Red']).add(1).divide(2
                            ).multiply(10000).rename('ndvi') #  escala [0, 1]

    ndwiImg = image.normalizedDifference(['Red', 'NIR']).add(1).divide(2
                            ).multiply(10000).rename('ndwi')   #  escala [0, 1]

    osaviImg = image.expression("float(b('NIR') - b('Red')) / (0.16 + b('NIR') + b('Red'))"
                    ).add(1).divide(2).multiply(10000).rename('osavi')  #  .add(1) escala [0, 1]  

    soilImg = image.expression("float(b('NIR') - b('Green')) / (b('NIR') + b('Green'))"
                    ).add(1).divide(2).multiply(10000).rename(['isoil'])    #  .add(1) escala [0, 1]
    
    # gcviImg = image.expression("float(b('NIR')) / (b('Green')) - 1").rename(['gcvi'])  
    
    # eviImg = eviImg.divide(10).add(1).divide(2).multiply(10000).rename('evi')  #escala [0, 1]  
    laiImg = laiImg.divide(10).add(1).divide(2).multiply(10000).rename('lai')  #escala [0, 1]  

    # ratioImg = image.expression("float(b('NIR') / b('Red'))").rename(['ratio'])  # escala [0, 1]

    image = image.addBands(laiImg).addBands(ndviImg)\
                .addBands(osaviImg).addBands(soilImg).addBands(ndwiImg)   # .addBands(eviImg)
                    
    # all bands in   escala [0, 1]
    return  image.select(bndInd).toUint16()


def agregateBandsgetFractions(image):

    outBandNames = ['gv', 'npv', 'soil', 'cloud','shade']
    
    # Define endmembers
    endmembers5 =  [
        [ 119.0,  475.0,  169.0, 6250.0, 2399.0,  675.0], #/*gv*/
        [1514.0, 1597.0, 1421.0, 3053.0, 7707.0, 1975.0], #/*npv*/
        [1799.0, 2479.0, 3158.0, 5437.0, 7707.0, 6646.0], #/*soil*/
        [4031.0, 8714.0, 7900.0, 8989.0, 7002.0, 6607.0], #/*cloud*/
        [   0.0,    0.0,    0.0,    0.0,    0.0,    0.0]  #/*Shade*/
    ]
    
    # Uminxing data  ==== Levando a a imagem a Intero primeiro de depois calculando 
    # as frações 
    fractions = ee.Image(image).multiply(10000).select([0,1,2,3,4,5])\
                    .unmix(endmembers= endmembers5, sumToOne= True, nonNegative= True).float()
    
    fractions = fractions.select([0,1,2,3,4], outBandNames).toFloat()

    image = image.addBands(fractions.select('soil')).addBands(
                        fractions.select('gv')).addBands(fractions.select('npv'))    
    return image.select(['soil', 'gv', 'npv'])

def agregateBandsIndexNDFIA(img):
    
    fracao = agregateBandsgetFractions(img)
    #calculate NDFI
    ndfia = fracao.expression(
        "float((b('gv') / (1 - b('shade'))) - b('soil')) / float((b('gv') / (1 - b('shade'))) + b('npv') + b('soil'))")
    
    ndfia = ndfia.add(1).multiply(255).toByte().rename('ndfia')    
    # return fracao.select(['soil','gv', 'npv']).addBands(ndfia).toFloat()
    return ndfia

def diferenciaImagens (imgA, imgB):    

    bndInd = ['lai','ndvi','ndwi','osavi','isoil'] 
    traslacao = [10000, 10000, 10000, 10000, 10000, 10000]
    bndVis = ['Blue','Green','Red','NIR']

    for cc, bnd in enumerate(bndInd):
        
        imgtemp = imgA.select(bnd).subtract(imgB.select(bnd)).add(
                                        traslacao[cc]).divide(2).toUint16().rename(bnd + '_d') 

        imgB = imgB.addBands(imgA.select([bnd], [bnd + '_m'])).addBands(imgtemp)

    for bnd in bndVis:
        imgB = imgB.addBands(imgA.select([bnd], [bnd + '_m']))
    
    return imgB.toUint16()


def image_correption_histogram_strecht (img):
    
    # refazendo os dados espectrais 
    img_HH = img.select(['Blue','Green','Red','NIR']).multiply(2.632).subtract(ee.Image.constant(526.4));
    
    # correção por baixo do zero 
    imgB2mask = img_HH.select('Blue').where(img_HH.select('Blue').lte(0), 0);
    imgB3mask = img_HH.select('Green').where(img_HH.select('Green').lte(0), 0);
    imgB4mask = img_HH.select('Red').where(img_HH.select('Red').lte(0), 0);
    imgB8mask = img_HH.select('NIR').where(img_HH.select('NIR').lte(0), 0);

    # correction por cima do zero 
    imgB2mask = imgB2mask.where(imgB2mask.gte(10000), 10000).rename('Blue');
    imgB3mask = imgB3mask.where(imgB3mask.gte(10000), 10000).rename('Green');
    imgB4mask = imgB4mask.where(imgB4mask.gte(10000), 10000).rename('Red');
    imgB8mask = imgB8mask.where(imgB8mask.gte(10000), 10000).rename('NIR');

    imgTransf = ee.Image.cat([imgB2mask,imgB3mask,imgB4mask,imgB8mask]).toUint16();

    return imgTransf


# ['blueBef', 'greenBef', 'redBef', 'NIRBef', 'blueAft', 'greenAft', 'redAft', 'NIRAft']
# Mean => [ 489.04966217  699.81953282  731.18343832 2380.71167809  504.76052985 706.04181667  735.5199312  2388.04778238]
# standard deviation => [329.22977418 366.72084416 443.58671308 662.45319986 361.94920024  395.08776032 467.65570749 702.11809735]
# Minimo  => [0, 0, 0, 0, 0, 0, 0, 0]
# Maximo => [20064.0, 18528.0, 17520.0, 16416.0, 23384.0, 20808.0, 19192.0, 16952.0]
#### ============================================================================ ####
# ['blueBef', 'greenBef', 'redBef', 'NIRBef', 'blueAft', 'greenAft', 'redAft', 'NIRAft']
# Mean => [ 398.11709736  628.99203297  673.56464596 2490.62616802  430.06746587  661.77794988  728.64307027 2544.52356961]
# standard deviation => [134.77702415 179.72047279 281.90808519 400.89204802 160.80124615 205.87923702 315.07080306 427.97224203]
# Minimo  => [0, 0, 0, 0, 0, 0, 0, 0]
# Maximo => [19312.0, 18208.0, 17184.0, 16320.0, 24600.0, 22008.0, 20024.0, 17896.0]


# https://www.kdnuggets.com/2020/04/data-transformation-standardization-normalization.html
def chips_MinMax_normalization_array(array_img, lst_min, lst_max):
    marray_img = copy.deepcopy(array_img)
    #(256,256, 8)
    for bnd in range(4):
        marray_img[:, :, bnd] = (array_img[:, :, bnd] - lst_min[bnd]) / (lst_max[bnd] - lst_min[bnd])
        marray_img[:, :, bnd + 4] = (array_img[:, :, bnd  + 4] - lst_min[bnd]) / (lst_max[bnd] - lst_min[bnd])
    
    return marray_img

# standardization (or Z-score normalization)
def chips_Standardisation_Zcore_array(array_img, lst_mean, lst_std):
    # xscaled = (x - means(x)) / std(x)
    marray_img = copy.deepcopy(array_img)
    #(256,256, 8)
    for bnd in range(4):
        marray_img[:, :, bnd] = (array_img[:, :, bnd] - lst_mean[bnd]) / lst_std[bnd] 
        marray_img[:, :, bnd + 4] = (array_img[:, :, bnd  + 4] - lst_mean[bnd]) / lst_std[bnd]
    
    return marray_img

def chips_Max_Min_scaling(array_img, lst_min, lst_max):
    # xscaled = (x - min(x)) / (max(x) - min(x))
    marray_img = copy.deepcopy(array_img)
    #(256,256, 8)
    for bnd in range(4):
        marray_img[:, :, bnd] = (array_img[:, :, bnd] - lst_min[bnd]) / (lst_max[bnd] - lst_min[bnd])
        marray_img[:, :, bnd + 4] = (array_img[:, :, bnd  + 4] - lst_min[bnd]) / (lst_max[bnd] - lst_min[bnd])
    
    return marray_img

def chips_Median_Quantiles_scaling(array_img, lst_median, lst_quantiles):
    # xscaled = (x - median(x)) / (75quantile(x) - 25quantile(x))
    marray_img = copy.deepcopy(array_img)
    #(256,256, 8)
    for bnd in range(4):
        marray_img[:, :, bnd] = (array_img[:, :, bnd] - lst_median[bnd]) / (lst_quantiles[bnd][1] - lst_quantiles[bnd][0])
        marray_img[:, :, bnd + 4] = (array_img[:, :, bnd + 4] - lst_median[bnd  + 4]) / (lst_quantiles[bnd + 4][1] - lst_quantiles[bnd  + 4][0])
    
    return marray_img

def calculing_zcore_modifing(arrX, xmean, xstd):
    calcZ = (arrX - xmean) / xstd
    return np.exp(-1 * calcZ)

def calculing_quartilesScore_modifing(arrX, lstQuart):
    calcZ = (arrX - lstQuart[1]) / (lstQuart[2] - lstQuart[0])
    return np.exp(-1 * calcZ)

def chips_Sigmoid_scaling(array_img, lst_mean, lst_std):
    # xscaled = 1 / ( 1 + exp(-(x- mean(x))/std(x)))
    marray_img = copy.deepcopy(array_img)
    #(256,256, 8)
    for bnd in range(4):
        expBandBef = calculing_zcore_modifing(array_img[:, :, bnd],lst_mean[bnd], lst_std[bnd])
        marray_img[:, :, bnd] = 1 / (1 + expBandBef)
        expBandAft = calculing_zcore_modifing(array_img[:, :, bnd  + 4],lst_mean[bnd  + 4], lst_std[bnd  + 4])
        marray_img[:, :, bnd + 4] =  1 / (1 + expBandAft)
    
    return marray_img

def chips_Robust_Sigmoid_scaling(array_img,  lst_quantiles):
    # xscaled = (x - min(x)) / (max(x) - min(x))
    marray_img = copy.deepcopy(array_img)
    #(256,256, 8)
    for bnd in range(4):
        expBandBef = calculing_quartilesScore_modifing(array_img[:, :, bnd], lst_quantiles[bnd])
        marray_img[:, :, bnd] = 1 / (1 + expBandBef)
        expBandAft = calculing_quartilesScore_modifing(array_img[:, :, bnd + 4], lst_quantiles[bnd  + 4])
        marray_img[:, :, bnd + 4] = 1 / (1 + expBandAft)
    
    return marray_img
