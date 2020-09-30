import ee


def apply(image, radius, units, time_window):

    time_window = dict(before=-time_window, after=time_window, units='month')
    t = image.date()
    fro = t.advance(ee.Number(time_window['before']), time_window['units'])
    to = t.advance(ee.Number(time_window['after']), time_window['units'])

    bands = image.bandNames()
    meanBands = bands.map(lambda b: ee.String(b).cat('_mean'))
    ratioBands = bands.map(lambda b: ee.String(b).cat('_ratio'))

    # create image collection for time-period of interest
    s1 = ee.ImageCollection('COPERNICUS/S1_GRD_FLOAT') \
        .filterMetadata('resolution_meters', 'equals', 10) \
        .filterMetadata('instrumentMode', 'equals', 'IW') \
        .filter(ee.Filter.And(
        ee.Filter.listContains('transmitterReceiverPolarisation', 'VV'),
        ee.Filter.listContains('transmitterReceiverPolarisation', 'VH')
        )) \
        .filterDate(fro, to) \
        .filter(
            ee.Filter.eq(
                'relativeOrbitNumber_start', image.get('relativeOrbitNumber_start')
            )
        )

    def mapMeanSpace(i):

        reducer = ee.Reducer.mean()
        kernel = ee.Kernel.square(radius, units)
        mean = i.reduceNeighborhood(reducer, kernel).rename(meanBands)
        ratio = i.divide(mean).rename(ratioBands)
        return i.addBands(mean).addBands(ratio).set()

    # compute spatial average for all images
    s1_mean = s1.map(mapMeanSpace)
    image_mean = s1_mean.filter(ee.Filter.eq('system:index', image.get('system:index'))).first()

    # get ratio bands
    s1_ratio = s1_mean.select(ratioBands)

    # do the filtering operation
    image_filtered = ee.Image(image_mean.select(meanBands)
        .multiply(s1_ratio.sum()).divide(s1_ratio.count())
        .rename(bands)
        .copyProperties(image)
        .set('system:time_start', image.get('system:time_start'))
    )

    return image_filtered
