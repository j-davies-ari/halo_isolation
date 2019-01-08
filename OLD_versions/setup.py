from distutils.core import setup, Extension

'''from Cython.Build import cythonizesetup(  name = 'Aperture creation app',  ext_modules = cythonize("apertures.pyx"),)
'''


module1 = Extension('apertures_cversion',                    sources = ['apertures_cversion.c'])setup (name = 'apertures_cversion',       version = '1.0',       description = 'yolo',       ext_modules = [module1])
