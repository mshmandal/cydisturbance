CycloneDisturbance
================
Mohammad Shamim Hasan Mandal
3/8/2021

The following script uses Landsat 5,7, and 8 NDVI images to compute
cyclone disturbances in the Sundarbans mangrove forest. The working
script in available here: [Google Earth
Engine](https://code.earthengine.google.com/4f2811375361dd9edad340ed22623889).
You need to have an GEE account to run the script.

Step 1: Import data layers to GEE

``` js
// Total sundarbans area
var aoi = ee.FeatureCollection("users/itzhasanshamim/Sundarban2018/Sundarban_Old_NatHaz");
// Classified image 
var sfc = ee.Image("users/itzhasanshamim/SFC");
//Simplified vegetation shapefile created using qgis
var veg = ee.FeatureCollection("users/itzhasanshamim/simplified_f_s_veg");                 

// Add classified vegetation map to map layer
Map.addLayer (veg,{}, 'Sundarban mangrove forest',1,0.3);

var ndviPalette = ['FFFFFF', 'CE7E45', 'DF923D', 'F1B555', 'FCD163', '99B718',
               '74A901', '66A000', '529400', '3E8601', '207401', '056201',
               '004C00', '023B01', '012E01', '011D01', '011301'];
               
Map.centerObject(aoi,9); //Center map to the area of interest

// Landsat 5, 7, 8 NDVI Composite
var l5col = ee.ImageCollection('LANDSAT/LT5_L1T_32DAY_NDVI');   // L5 NDVI
var l7col = ee.ImageCollection('LANDSAT/LE7_L1T_32DAY_NDVI');   // L7 NDVI
var l8col = ee.ImageCollection('LANDSAT/LC8_L1T_32DAY_NDVI');   // L8 NDVI
```

Step 2: Make data variable for filtering the image collection

Images between Jan-March will be used for this analysis. However, for
the cyclone event of 2006 and 2008, this range is extened upto May to
produce relevant image mosaic.

``` js
var y88 = ee.Filter.date('1988-01-01','1988-03-31');
var y89 = ee.Filter.date('1989-01-01','1989-03-31'); 
var y90 = ee.Filter.date('1990-01-01','1990-03-31');
var y91 = ee.Filter.date('1991-01-01','1991-03-31');
var y92 = ee.Filter.date('1992-01-01','1992-03-31');
var y93 = ee.Filter.date('1993-01-01','1993-03-31'); // No cylone
var y94 = ee.Filter.date('1994-01-01','1994-03-31');
var y95 = ee.Filter.date('1995-01-01','1995-03-31');
var y96 = ee.Filter.date('1996-01-01','1996-03-31');
var y97 = ee.Filter.date('1997-01-01','1997-03-31');
var y98 = ee.Filter.date('1998-01-01','1998-03-31');
var y99 = ee.Filter.date('1999-01-01','1999-03-31');
var y00 = ee.Filter.date('2000-01-01','2000-03-31');
var y01 = ee.Filter.date('2001-01-01','2001-03-31'); // No cylone
var y02 = ee.Filter.date('2002-01-01','2002-03-31'); 
var y03 = ee.Filter.date('2003-01-01','2003-03-31'); // No cylone
var y04 = ee.Filter.date('2004-01-01','2004-03-31');
var y05 = ee.Filter.date('2005-01-01','2005-03-31');
var y06 = ee.Filter.date('2006-01-01','2006-03-31'); 
var y07s = ee.Filter.date('2007-01-01','2007-04-28'); // For Post C06
var y07 = ee.Filter.date('2007-01-01','2007-03-31');
var y08 = ee.Filter.date('2008-01-01','2008-03-31');
var y08s = ee.Filter.date('2008-01-01','2008-05-31'); // For Pre C08;
var y09 = ee.Filter.date('2009-01-01','2009-03-31');
var y10 = ee.Filter.date('2010-01-01','2010-03-31'); // No cylone
var y15 = ee.Filter.date('2015-01-01','2015-03-31');
var y16 = ee.Filter.date('2016-01-01','2016-03-31');
var y17 = ee.Filter.date('2017-01-01','2017-03-31'); // No cylone
```

Step 3: Create functions

As we want to compute forest area disturbances by multiple cyclones. We
will create few functions to apply them in each events.

Function 1: Filter function \[Make functions to filter and mask\]

``` js
var filterFun = function (imageCollection, dateRange){
  var fCol = imageCollection.filter(dateRange)
              // Get year and add it to the image property
              .map(function(img) {
              return img.set('year', ee.Image(img).date().get('year'));
              });
  
  // Appply quality mosaic using NDVI and clip AOI
  var newImage= fCol.qualityMosaic('NDVI').clip(aoi)
                    .set('year', fCol.first().get('year'));

  var waterMask = newImage.gte(0);
  var filteredImage = newImage.mask(waterMask);
  //Map.addLayer (filteredImage, ndviPalette, "filteredImage mask",0);
  return filteredImage;
};
```

Function 2: function to computer pre and post event disturbances \[using
thersholding\]

``` js
var lossFun = function (preImage, postImage){
  var difference = preImage.subtract(postImage);   // Subtract preImage from postImage
  var diminished= difference.gt(0.1);              // NDVI changes more than 0.1
  var masked = diminished.mask(diminished);        // Mask the areas that meet the threshold condition
  var mask = masked.clip(aoi);                     // Clip to area of interest
  // Masking out non forest area
  var fmask = sfc.select('classification').eq(0);
  var ndviMasked = preImage.updateMask(fmask);
  //Map.addLayer (ndviMasked, ndviPalette, "Forest mask",0); 
  var area = ee.Image.pixelArea().divide(1000000); // variable for converting "sq meters" to "square kilometers";
  var preNDVI = preImage.multiply(area);           // Forest area, before the event
  var ndviAffectedArea = mask.multiply(area);      // Forest area, disturbed
  var areas = ndviAffectedArea.rename('Lost_NDVI') // Save the disturbed area and rename
            .addBands(preNDVI.rename('Pre_NDVI')); // Save the area before event and rename
            
  return areas;
};
```

Function 3: Calculate disturbed area, and area before disturbance

``` js
var areaFun = function(image){
  var x = image.reduceRegion(
      {
          'reducer': ee.Reducer.sum(),
          'scale':30,
          'geometry': veg,  // Use the geometry [vegetation area, from classified image] to calculate 
          'maxPixels': 1e12
        });
        return x;
};
```

Funtion 4: Calculate for each events

``` js
var eventFun = function(imageCol, preYear, postYear){
  var pre = filterFun(imageCol, preYear); 
  var post = filterFun(imageCol, postYear);
  // Calculate area in pre event
  var area = ee.Image.pixelArea().divide(1000000);
  var preCal = pre.multiply(area);
  
  var preArea_f = areaFun(preCal);              // Pre event area
  //print("Area_before_disturbance:", preArea_f); 
  
  var m_pre = pre.updateMask(sfc.select("classification").eq(0));
  var m_post = post.updateMask(sfc.select("classification").eq(0));
  var eventYear = ee.String(pre.get('year')).getInfo();
  //print('eventYear: ', eventYear); // ee.Date
  
  //Map.addLayer(m_pre, {'palette': ndviPalette, 'min': 0, 'max': 1}, 'Pre event',0);
  //Map.addLayer(m_post,{'palette': ndviPalette, 'min': 0, 'max': 1}, 'Post event',0);
  var lossBand_f = lossFun(pre, post);            // apply "lossFun" function
  //print("Loss band from function_F4",lossBand_f);  // apply "lossFun" function
  Map.addLayer(lossBand_f.select("Lost_NDVI"), {palette: ['FF0000'],min: 0, max: 1}, 'Loss in the year '+eventYear,0); 
  // Estimate lost area
  var lostArea_f = areaFun(lossBand_f);          // apply "areaFun" function
```

Error check : START

This code chuck was used to recheck the area calculation.

``` js
// Estimation of pre and post area 
var pixelArea = ee.Image.pixelArea().divide(1000000);
var preEvent_unmasked = pre.multiply(pixelArea);
var postEvent_unmasked = post.multiply(pixelArea);

var preEvent_image = m_pre.multiply(pixelArea);
var postEvent_image = m_post.multiply(pixelArea);

  
var getAreas = preEvent_unmasked.rename("1_pre_unmasked").addBands(postEvent_unmasked.rename("2_post_unmasked")).addBands(preEvent_image
                            .rename("3_pre_masked")).addBands(postEvent_image.rename("4_post_masked"));
                            
var ErrorCheck_Stats = getAreas.reduceRegion(
  {       'reducer': ee.Reducer.sum(),
          'scale':30,
          'geometry': aoi,
          'maxPixels': 1e9
        });

//print("ErrorCheck_Stats:", ErrorCheck_Stats);
```

Error check : END

``` js
  // Estimate total forest area from the classified image 
  var totalF1 = sfc.select("classification").eq(0);
  var totalF2 = totalF1.multiply(area);
  var totalF3 = areaFun(totalF2);               // apply "areaFun" function
  
  
  //print("Loss_are and Pre NDVI area_F4",lostArea_f);
  //print("Area Before_F4",preArea_f);
  //print("Total forest from Classification?,",totalF3);
  
// Extra -----check area results

  // return values
 return {
        preImage: pre,
        postImage: post,
        lossBand: lossBand_f,
        lossArea: lostArea_f,
    }; 
  };
```

Step 4: Apply the functions to analyze the filtered image collection

Create separate variables from function outputs.

``` js
// C88
var c88 = eventFun(l5col,y88, y89); // Landsat 5

// C89
var c89 = eventFun(l5col,y89, y90); // Landsat 5

// C90
var c90 = eventFun(l5col,y90, y91); // Landsat 5

// C91
var c91 = eventFun(l5col,y91, y92); // Landsat 5

// C92
var c92 = eventFun(l5col,y92, y93); // Landsat 5

// C94
var c94 = eventFun(l5col,y94, y95); // Landsat 5

// C95
var c95 = eventFun(l5col,y95, y96); // Landsat 5

// C96
var c96 = eventFun(l5col,y96, y97); // Landsat 5

// C97
var c97 = eventFun(l5col,y97, y98); // Landsat 5

// C98
var c98 = eventFun(l5col,y98, y99); // Landsat 5

// C99
var c99 = eventFun(l5col,y99, y00); // Landsat 5

// C00
var c00 = eventFun(l5col,y00, y01); // Landsat 5

// C02
var c02 = eventFun(l7col,y02, y03); // Landsat 7

//C04
var c04 = eventFun(l5col,y04, y05); // Landsat 5

// C05
var c05 = eventFun(l5col,y05, y06); // Landsat 5

// C06
var c06 = eventFun(l5col,y06, y07s); // Landsat 5

// C07
var c07 = eventFun(l5col,y07, y08); // Landsat 5

// C08
var c08 = eventFun(l5col,y08s, y09); // Landsat 5

// C09
var c09 = eventFun(l5col,y09, y10); // Landsat 5

// C015
var c15 = eventFun(l8col,y15, y16); // Landsat 8

// C16
var c16 = eventFun(l8col,y16, y17); // Landsat 8
```

Make a Dictionary of lost NDVI per cyclone

``` js
var dict = ee.Dictionary({
c88: c88.lossArea.get('Lost_NDVI'),
c89: c89.lossArea.get('Lost_NDVI'),
c90: c90.lossArea.get('Lost_NDVI'), 
c91: c91.lossArea.get('Lost_NDVI'), 
c92: c92.lossArea.get('Lost_NDVI'), 
c94: c94.lossArea.get('Lost_NDVI'),
c95: c95.lossArea.get('Lost_NDVI'), 
c96: c96.lossArea.get('Lost_NDVI'), 
c97: c97.lossArea.get('Lost_NDVI'),
c98: c98.lossArea.get('Lost_NDVI'), 
c99: c99.lossArea.get('Lost_NDVI'),
c00: c00.lossArea.get('Lost_NDVI'),
c02: c02.lossArea.get('Lost_NDVI'), 
c04: c04.lossArea.get('Lost_NDVI'),
c05: c05.lossArea.get('Lost_NDVI'),
c06: c06.lossArea.get('Lost_NDVI'),
c07: c07.lossArea.get('Lost_NDVI'), 
c08: c08.lossArea.get('Lost_NDVI'), 
c09: c09.lossArea.get('Lost_NDVI'), 
c15: c15.lossArea.get('Lost_NDVI'), 
c16: c16.lossArea.get('Lost_NDVI')
});

// Print the dictionary
print("Affected forest area (sq km), in terms of NDVI loss: ",dict);
```

End of script

Project: Mandal, M.S.H., Hosaka, T. Assessing cyclone disturbances
(1988–2016) in the Sundarbans mangrove forests using Landsat and Google
Earth Engine. Nat Hazards 102, 133–150 (2020).
<https://doi.org/10.1007/s11069-020-03914-z>
