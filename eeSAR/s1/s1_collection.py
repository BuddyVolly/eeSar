import ee
import math
from functools import reduce

from ..misc import dates
from ..misc import _refined_lee
from ..misc import _quegan
from ..misc import _slope_correction


def create(
        region,
        orbits=['ASCENDING', 'DESCENDING'],
        start_date=None,
        end_date=None,
        target_date=None,
        add_ratio=True,
        add_ND_ratio=True,
        speckle_filter='NONE',
        speckle_filter_dict={'radius': 30, 'units': 'meters',
                             'time_window': 12},
        radiometric_correction='ELLIPSOID',
        slope_correction_dict={'model': 'volume', 'dem': 'USGS/SRTMGL1_003',
                               'buffer': 12},
        outlier_removal='NONE',
        db=True
):
    """Creates an ee.ImageCollection with Sentinel 1 imagery.

    Args:
        region: The region to include imagery within. (ee.Geometry)

        orbits:  The orbits to include. (string: ASCENDING and/or DESCENDING)

        start_date: (Optional) The earliest date to include images for (inclusive).
            Defaults to 6 months before target date if target date is specified,
            Can be a string (yyyy-mm-dd), python datetime, or millis timestamp.

        end_date: (Optional) The latest date to include images for (exclusive).
            Defaults to 6 months after target date if target date is specified,
            Can be a String (yyyy-mm-dd), python datetime, or millis timestamp.

        target_date: (Optional) The ideal acquisition date of imagery.

        speckle_filter: (Optional) Type of speckle filtering to apply. Can be one of the following:
            NONE - No speckle filtering is applied;
            BOXCAR - 30x30m boxcar filtering is applied;
            SNIC - Simple Non-Iterative Clustering is used to filter out speckle;
            REFINED_LEE - Refined Lee filtering is applied. This can be very slow, but often gives good results.
            QUEGAN - Multitemporal Speckle filter. Slow but good.

        speckle_filter_dict: (Optional) Parameters for Quegan's MT Speckle Filter
            radius -  radius for the mean filter
            units - units for the radius
            time_window (in months) - images to consider before and after that time-window

            parameters as given as dict, e.g.:
            {'radius': 30, 'units': 'meters', 'time_window': 12}

        radiometric_correction: (Optional) The geometric correction to apply. Can be one of the following:
            NONE - No correction is applied;
            ELLIPSOID (default) - Ellipsoid correction (gamma0 correction) is applied;
            TERRAIN - Terrain correction is applied.

        slope_correction_dict: (Optional) Parameter for radiometric TERRAIN Correction
            model - which model to use (surface, or volume)
            dem - which elevation model to use (as EE asset)
            buffer - buffer parameter for layover/shadow mask

        outlier_removal: (Optional) Removes outliers from the image collection. Can be one of the following:
            NONE - No outlier removal is performed;
            MODERATE - Remove outliers
            AGGRESSIVE - Remove more outliers

    Returns:
        an ee.ImageCollection with Sentinel 1 imagery with the following band:
            VV;
            VH;
            angle.

        If targetDate is specified, the following metadata bands are also included:
            dayOfYear - day of the year;
            daysFromTarget - days from the specified target date;
            quality - closeness to target date;
            unixTimeDays - days since unix epoch.
    """

    days = 366 / 2
    if target_date and not start_date:
        start_date = dates.subtract_days(target_date, days)

    if target_date and not end_date:
        end_date = dates.add_days(target_date, days)

    def set_resample(image):
        """ Set resampling of the image to bilinear"""
        return image.resample()

    def to_gamma0(image):
        """ Calibrate to ellipsoid corrected Gamma0

        This bases on SAR image data in power scale

        :param image:
        :return:
        """
        gamma0 = image.expression('i/(cos(angle * pi / 180))', {
            'i': image.select(['VV', 'VH']),
            'angle': image.select('angle'),
            'pi': math.pi
        })
        return image.addBands(gamma0, None, True)

    def boxcar_filter(image):

        filtered = image.select(['VV', 'VH']) \
            .reduceNeighborhood(ee.Reducer.mean(),
                                ee.Kernel.square(30, 'meters')) \
            .rename(['VV', 'VH'])
        return image.addBands(filtered, None, True)

    def quegan(image):

        radius = speckle_filter_dict['radius']
        units = speckle_filter_dict['units']
        time_window = speckle_filter_dict['time_window']

        return _quegan.apply(image, radius, units, time_window)

    def terrain_flattening(image):

        model = slope_correction_dict['model']
        dem = slope_correction_dict['dem']
        buffer = slope_correction_dict['buffer']

        return _slope_correction.correct(image, model, dem, buffer)

    def mask_overlay(image):
        return image.updateMask(
            image.select('layover').And(image.select('shadow'))
        )

    def snic_filter(image):

        bands = ['VV', 'VH', 'VVVH_ratio', 'VVVH_ND_ratio']
        snic = ee.Algorithms.Image.Segmentation.SNIC(
            image=image.select(bands),
            size=8,
            compactness=5,
            connectivity=8,
            neighborhoodSize=16,
        ).select('.*_mean').rename(bands)

        return image.addBands(snic, None, True)

    def add_date_bands(image):
        date = image.date()
        day_of_year = date.getRelative('day', 'year')
        days_from_target = date.difference(ee.Date(target_date), 'day').abs()
        millisPerDay = 1000 * 60 * 60 * 24
        unix_time_days = date.millis().divide(millisPerDay)
        return image \
            .addBands(ee.Image(day_of_year).uint16().rename('dayOfYear')) \
            .addBands(
                ee.Image(days_from_target).uint16().rename('daysFromTarget')) \
            .addBands(
                ee.Image(days_from_target.multiply(-1)).int16().rename('quality')) \
            .addBands(ee.Image(unix_time_days).uint16().rename('unixTimeDays'))

    def add_ratio_band(image):
        return image \
            .addBands(image.select('VV').subtract(image.select('VH')) \
                      .rename('VVVH_ratio'))

    def add_ND_ratio_band(image):
        return image \
            .addBands(image.normalizedDifference(['VV', 'VH']) \
                      .rename('VVVH_ND_ratio'))

    def mask_borders(image):
        total_slices = ee.Number(image.get('totalSlices'))
        slice_number = ee.Number(image.get('sliceNumber'))
        middle_slice = ee.Image(
            slice_number.gt(1).And(slice_number.lt(total_slices)))
        not_border = image.select('VV').mask().reduceNeighborhood(
            reducer=ee.Reducer.allNonZero(),
            kernel=ee.Kernel.circle(radius=500, units='meters')
        )
        angle = image.select('angle')
        return image \
            .updateMask(
            angle.gt(31).And(angle.lt(45)).And(middle_slice.Or(not_border))
        )

    def toDb(image):

        db = ee.Image(10).multiply(image.select(['VV', 'VH']).log10()) \
            .rename(['VV', 'VH'])
        return image.addBands(db, None, True)

    def pre_process(image):
        steps = []

        if speckle_filter == 'BOXCAR':
            steps.append(boxcar_filter)
        elif speckle_filter == 'REFINED_LEE':
            steps.append(_refined_lee.apply)
        elif speckle_filter == 'QUEGAN':
            steps.append(quegan)
        if radiometric_correction == 'ELLIPSOID':
            steps.append(to_gamma0)
        elif radiometric_correction == 'TERRAIN':
            steps.append(terrain_flattening)
            steps.append(mask_overlay)
        if target_date:
            steps.append(add_date_bands)
        if add_ratio:
            steps.append(add_ratio_band)
        if add_ND_ratio:
            steps.append(add_ND_ratio_band)
        if speckle_filter == 'SNC':
            steps.append(snic_filter)
        if db:
            steps.append(toDb)

        steps.append(mask_borders)



        return reduce(lambda acc, process: process(acc), steps, image)

    def mask_outliers(collection, std_devs):
        bands = ['VV', 'VH']
        reduced = collection.select(bands).reduce(
            ee.Reducer.median().combine(ee.Reducer.stdDev(), '', True)
        )
        median = reduced.select('.*_median')
        std_dev = reduced.select('.*_stdDev')
        threshold = std_dev.multiply(std_devs)
        maskedCollection = collection.map(
            lambda image: image.updateMask(
                image.select(bands).subtract(median).abs().lte(threshold)
                    .reduce(ee.Reducer.min())
            )
        )
        return maskedCollection

    collection = ee.ImageCollection('COPERNICUS/S1_GRD_FLOAT') \
        .filterBounds(region) \
        .filterDate(start_date, end_date) \
        .filterMetadata('resolution_meters', 'equals', 10) \
        .filter(ee.Filter.eq('instrumentMode', 'IW')) \
        .filter(ee.Filter.And(
        ee.Filter.listContains('transmitterReceiverPolarisation', 'VV'),
        ee.Filter.listContains('transmitterReceiverPolarisation', 'VH')
    )) \
        .filter(ee.Filter.Or(
        [ee.Filter.eq('orbitProperties_pass', orbit) for orbit in orbits]
    )) \
        .map(set_resample)\
        .map(pre_process)\

    if outlier_removal == 'MODERATE':
        return mask_outliers(collection, 3)
    elif outlier_removal == 'AGGRESSIVE':
        return mask_outliers(collection, 2.6)
    else:
        return collection
