# CEBALERT: incomplete - being written.

# To do:
# - finish the positioning tests (like tests 1 and 2)
# - test when image is larger than sheet
# - test edge_average
# - test rotation
# - test repositioning


### You might want to use topo.command.pylabplot.matrixplot to
### visualize matrices during debugging.

import unittest
from numpy.testing import assert_array_almost_equal

from numpy.oldnumeric import array,Float,pi
## from numpy.oldnumeric.mlab import rot90

from param import resolve_path

from topo.base.boundingregion import BoundingBox
from topo.pattern.image import FileImage
from topo.transferfn import IdentityTF



class TestImage(unittest.TestCase):

    def test_vertical_oddimage_evensheet__horizontal_evenimage_evensheet(self):
        """
        Test vertical positioning for even sheet, odd image and horizontal
        positioning for even image, even sheet.
        """
        image_array = array(
[[  96.59090909,  96.59090909,  96.59090909,  96.59090909,  96.59090909,
         96.59090909,  96.59090909,  96.59090909,],
 [   0.        ,  34.        ,  68.        , 102.        , 136.        ,
        255.   ,   0.        ,   0.             ,],
 [   0.        ,  34.        ,  68.        , 102.        , 136.        ,
        255.        , 255.        ,   0.        ,],
 [   0.        ,  34.        ,  68.        , 102.        , 136.        ,
        255.        , 255.        , 255.        ,],
 [ 255.        ,   0.        , 255.        ,   0.        , 255.        ,
          0.        , 255.        ,   0.        ,],
 [   0.        , 255.        ,   0.        , 255.        ,   0.        ,
        255.        ,   0.        , 255.        ,],
 [  96.59090909,  96.59090909,  96.59090909,  96.59090909,  96.59090909,
         96.59090909,  96.59090909,  96.59090909,],
 [  96.59090909,  96.59090909,  96.59090909,  96.59090909,  96.59090909,
         96.59090909,  96.59090909,  96.59090909,]],Float)


        image = FileImage(filename = resolve_path('topo/tests/unit/testimage.pgm'),
                      xdensity=8,
                      ydensity=8,
                      bounds=BoundingBox(radius=0.5),
                      output_fns=[])

        ps = image.pattern_sampler
        ps.size_normalization='original'
        ps.whole_pattern_output_fns=[]

        assert_array_almost_equal(image_array,image())



    # CEBHACKALERT: this test is failing because of sheet2matrix()/sheet2matrixidx();
    # it might not be possible to stop it failing. It's skipped for the moment.
    def test_vertical_oddimage_oddsheet__horizontal_evenimage_oddsheet(self):
        """
        Test vertical centering for odd sheet, odd image, and horizontal
        centering for odd sheet, even image.

        FileImage is smaller than Sheet on which it's displayed.
        """
        pass
##         image_array = array(
## [[  96.59090909,  96.59090909,  96.59090909,  96.59090909,  96.59090909,
##          96.59090909,  96.59090909,  96.59090909,  96.59090909,],
##  [  96.59090909,  96.59090909,  96.59090909,  96.59090909,  96.59090909,
##          96.59090909,  96.59090909,  96.59090909,  96.59090909,],
##  [   0.        ,  34.        ,  68.        , 102.        , 136.        ,
##         255.   ,   0.        ,   0.             ,  96.59090909,],
##  [   0.        ,  34.        ,  68.        , 102.        , 136.        ,
##         255.        , 255.        ,   0.        ,  96.59090909,],
##  [   0.        ,  34.        ,  68.        , 102.        , 136.        ,
##         255.        , 255.        , 255.        ,  96.59090909,],
##  [ 255.        ,   0.        , 255.        ,   0.        , 255.        ,
##           0.        , 255.        ,   0.        ,  96.59090909,],
##  [   0.        , 255.        ,   0.        , 255.        ,   0.        ,
##         255.        ,   0.        , 255.        ,  96.59090909,],
##  [  96.59090909,  96.59090909,  96.59090909,  96.59090909,  96.59090909,
##          96.59090909,  96.59090909,  96.59090909,  96.59090909,],
##  [  96.59090909,  96.59090909,  96.59090909,  96.59090909,  96.59090909,
##          96.59090909,  96.59090909,  96.59090909,  96.59090909,]],Float)


##         image = FileImage(filename = resolve_path('topo/tests/unit/testimage.pgm'),
##                       xdensity=9,
##                       ydensity=9,
##                       output_fn=IdentityTF(),
##                       whole_image_output_fn=IdentityTF(),
##                       size_normalization='original')

##         assert_array_almost_equal(image_array,image())



    def test_fit_longest(self):
        """
        Check that the longer image dimension is made to fit the default
        dimension of 1.0, while the other is scaled the same.
        """
        ### Twice the default BoundingBox dimensions, image size of 2.0.
        ### In this case, 8 units represent 1.0 in sheet coordinates.
        image_array = array(
[[  96.59090909,  96.59090909,  96.59090909,  96.59090909,  96.59090909,
         96.59090909,  96.59090909,  96.59090909,  96.59090909,  96.59090909,
         96.59090909,  96.59090909,  96.59090909,  96.59090909,  96.59090909,
         96.59090909,],
 [  96.59090909,  96.59090909,  96.59090909,  96.59090909,  96.59090909,
         96.59090909,  96.59090909,  96.59090909,  96.59090909,  96.59090909,
         96.59090909,  96.59090909,  96.59090909,  96.59090909,  96.59090909,
         96.59090909,],
 [  96.59090909,  96.59090909,  96.59090909,  96.59090909,  96.59090909,
         96.59090909,  96.59090909,  96.59090909,  96.59090909,  96.59090909,
         96.59090909,  96.59090909,  96.59090909,  96.59090909,  96.59090909,
         96.59090909,],
 [   0.        ,   0.        ,  34.        ,  34.        ,  68.        ,
         68.        , 102.        , 102.        , 136.        , 136.        ,
        255.        , 255.        ,   0.        ,   0.        ,   0.        ,
          0.        ,],
 [   0.        ,   0.        ,  34.        ,  34.        ,  68.        ,
         68.        , 102.        , 102.        , 136.        , 136.        ,
        255.        , 255.        ,   0.        ,   0.        ,   0.        ,
          0.        ,],
 [   0.        ,   0.        ,  34.        ,  34.        ,  68.        ,
         68.        , 102.        , 102.        , 136.        , 136.        ,
        255.        , 255.        , 255.        , 255.        ,   0.        ,
          0.        ,],
 [   0.        ,   0.        ,  34.        ,  34.        ,  68.        ,
         68.        , 102.        , 102.        , 136.        , 136.        ,
        255.        , 255.        , 255.        , 255.        ,   0.        ,
          0.        ,],
 [   0.        ,   0.        ,  34.        ,  34.        ,  68.        ,
         68.        , 102.        , 102.        , 136.        , 136.        ,
        255.        , 255.        , 255.        , 255.        , 255.        ,
        255.        ,],
 [   0.        ,   0.        ,  34.        ,  34.        ,  68.        ,
         68.        , 102.        , 102.        , 136.        , 136.        ,
        255.        , 255.        , 255.        , 255.        , 255.        ,
        255.        ,],
 [ 255.        , 255.        ,   0.        ,   0.        , 255.        ,
        255.        ,   0.        ,   0.        , 255.        , 255.        ,
          0.        ,   0.        , 255.        , 255.        ,   0.        ,
          0.        ,],
 [ 255.        , 255.        ,   0.        ,   0.        , 255.        ,
        255.        ,   0.        ,   0.        , 255.        , 255.        ,
          0.        ,   0.        , 255.        , 255.        ,   0.        ,
          0.        ,],
 [   0.        ,   0.        , 255.        , 255.        ,   0.        ,
          0.        , 255.        , 255.        ,   0.        ,   0.        ,
        255.        , 255.        ,   0.        ,   0.        , 255.        ,
        255.        ,],
 [   0.        ,   0.        , 255.        , 255.        ,   0.        ,
          0.        , 255.        , 255.        ,   0.        ,   0.        ,
        255.        , 255.        ,   0.        ,   0.        , 255.        ,
        255.        ,],
 [  96.59090909,  96.59090909,  96.59090909,  96.59090909,  96.59090909,
         96.59090909,  96.59090909,  96.59090909,  96.59090909,  96.59090909,
         96.59090909,  96.59090909,  96.59090909,  96.59090909,  96.59090909,
         96.59090909,],
 [  96.59090909,  96.59090909,  96.59090909,  96.59090909,  96.59090909,
         96.59090909,  96.59090909,  96.59090909,  96.59090909,  96.59090909,
         96.59090909,  96.59090909,  96.59090909,  96.59090909,  96.59090909,
         96.59090909,],
 [  96.59090909,  96.59090909,  96.59090909,  96.59090909,  96.59090909,
         96.59090909,  96.59090909,  96.59090909,  96.59090909,  96.59090909,
         96.59090909,  96.59090909,  96.59090909,  96.59090909,  96.59090909,
         96.59090909,]],Float)


        image = FileImage(filename = resolve_path('topo/tests/unit/testimage.pgm'),
                      xdensity=8,
                      ydensity=8,
                      size=2.0,
                      output_fns=[],
                      bounds=BoundingBox(radius=1.0))

        ps = image.pattern_sampler
        ps.size_normalization='fit_longest'
        ps.whole_pattern_output_fns=[]

        assert_array_almost_equal(image_array,image())



    def test_stretch_to_fit(self):
        """
        Test that both image dimensions are made to fit 1.0.
        """
        ### 8 units represent 1.0 in sheet coordinates.
        image_array = array(
[[  96.59090909,  96.59090909,  96.59090909,  96.59090909,  96.59090909,
         96.59090909,  96.59090909,  96.59090909,  96.59090909,  96.59090909,
         96.59090909,  96.59090909,  96.59090909,  96.59090909,  96.59090909,
         96.59090909,],
 [  96.59090909,  96.59090909,  96.59090909,  96.59090909,  96.59090909,
         96.59090909,  96.59090909,  96.59090909,  96.59090909,  96.59090909,
         96.59090909,  96.59090909,  96.59090909,  96.59090909,  96.59090909,
         96.59090909,],
 [  96.59090909,  96.59090909,  96.59090909,  96.59090909,  96.59090909,
         96.59090909,  96.59090909,  96.59090909,  96.59090909,  96.59090909,
         96.59090909,  96.59090909,  96.59090909,  96.59090909,  96.59090909,
         96.59090909,],
 [  96.59090909,  96.59090909,  96.59090909,  96.59090909,  96.59090909,
         96.59090909,  96.59090909,  96.59090909,  96.59090909,  96.59090909,
         96.59090909,  96.59090909,  96.59090909,  96.59090909,  96.59090909,
         96.59090909,],
 [  96.59090909,  96.59090909,  96.59090909,  96.59090909,   0.        ,
         34.        ,  68.        , 102.        , 136.        , 255.        ,
          0.        ,   0.        ,  96.59090909,  96.59090909,  96.59090909,
         96.59090909,],
 [  96.59090909,  96.59090909,  96.59090909,  96.59090909,   0.        ,
         34.        ,  68.        , 102.        , 136.        , 255.        ,
          0.        ,   0.        ,  96.59090909,  96.59090909,  96.59090909,
         96.59090909,],
 [  96.59090909,  96.59090909,  96.59090909,  96.59090909,   0.        ,
         34.        ,  68.        , 102.        , 136.        , 255.        ,
        255.        ,   0.        ,  96.59090909,  96.59090909,  96.59090909,
         96.59090909,],
 [  96.59090909,  96.59090909,  96.59090909,  96.59090909,   0.        ,
         34.        ,  68.        , 102.        , 136.        , 255.        ,
        255.        , 255.        ,  96.59090909,  96.59090909,  96.59090909,
         96.59090909,],
 [  96.59090909,  96.59090909,  96.59090909,  96.59090909,   0.        ,
         34.        ,  68.        , 102.        , 136.        , 255.        ,
        255.        , 255.        ,  96.59090909,  96.59090909,  96.59090909,
         96.59090909,],
 [  96.59090909,  96.59090909,  96.59090909,  96.59090909, 255.        ,
          0.        , 255.        ,   0.        , 255.        ,   0.        ,
        255.        ,   0.        ,  96.59090909,  96.59090909,  96.59090909,
         96.59090909,],
 [  96.59090909,  96.59090909,  96.59090909,  96.59090909,   0.        ,
        255.        ,   0.        , 255.        ,   0.        , 255.        ,
          0.        , 255.        ,  96.59090909,  96.59090909,  96.59090909,
         96.59090909,],
 [  96.59090909,  96.59090909,  96.59090909,  96.59090909,   0.        ,
        255.        ,   0.        , 255.        ,   0.        , 255.        ,
          0.        , 255.        ,  96.59090909,  96.59090909,  96.59090909,
         96.59090909,],
 [  96.59090909,  96.59090909,  96.59090909,  96.59090909,  96.59090909,
         96.59090909,  96.59090909,  96.59090909,  96.59090909,  96.59090909,
         96.59090909,  96.59090909,  96.59090909,  96.59090909,  96.59090909,
         96.59090909,],
 [  96.59090909,  96.59090909,  96.59090909,  96.59090909,  96.59090909,
         96.59090909,  96.59090909,  96.59090909,  96.59090909,  96.59090909,
         96.59090909,  96.59090909,  96.59090909,  96.59090909,  96.59090909,
         96.59090909,],
 [  96.59090909,  96.59090909,  96.59090909,  96.59090909,  96.59090909,
         96.59090909,  96.59090909,  96.59090909,  96.59090909,  96.59090909,
         96.59090909,  96.59090909,  96.59090909,  96.59090909,  96.59090909,
         96.59090909,],
 [  96.59090909,  96.59090909,  96.59090909,  96.59090909,  96.59090909,
         96.59090909,  96.59090909,  96.59090909,  96.59090909,  96.59090909,
         96.59090909,  96.59090909,  96.59090909,  96.59090909,  96.59090909,
         96.59090909,]])


        image = FileImage(filename = resolve_path('topo/tests/unit/testimage.pgm'),
                      xdensity=8,
                      ydensity=8,
                      output_fns=[],
                      bounds=BoundingBox(radius=1.0))

        ps = image.pattern_sampler
        ps.size_normalization='stretch_to_fit'
        ps.whole_pattern_output_fns=[]

        assert_array_almost_equal(image_array,image())



    def test_fit_shortest(self):
        """
        Test that the shorter dimension is made to fit 1.0, while the other
        is scaled by the same factor.
        """
        ### 15 units represent 1.0 in sheet coordinates.
        image_array = array(
 [[  34.,  68.,  68.,  68., 102., 102., 102., 136., 136., 136., 255., 255., 255.,
            0.,   0.,],
 [  34.,  68.,  68.,  68., 102., 102., 102., 136., 136., 136., 255., 255., 255.,
            0.,   0.,],
 [  34.,  68.,  68.,  68., 102., 102., 102., 136., 136., 136., 255., 255., 255.,
            0.,   0.,],
 [  34.,  68.,  68.,  68., 102., 102., 102., 136., 136., 136., 255., 255., 255.,
          255., 255.,],
 [  34.,  68.,  68.,  68., 102., 102., 102., 136., 136., 136., 255., 255., 255.,
          255., 255.,],
 [  34.,  68.,  68.,  68., 102., 102., 102., 136., 136., 136., 255., 255., 255.,
          255., 255.,],
 [  34.,  68.,  68.,  68., 102., 102., 102., 136., 136., 136., 255., 255., 255.,
          255., 255.,],
 [  34.,  68.,  68.,  68., 102., 102., 102., 136., 136., 136., 255., 255., 255.,
          255., 255.,],
 [  34.,  68.,  68.,  68., 102., 102., 102., 136., 136., 136., 255., 255., 255.,
          255., 255.,],
 [   0., 255., 255., 255.,   0.,   0.,   0., 255., 255., 255.,   0.,   0.,   0.,
          255., 255.,],
 [   0., 255., 255., 255.,   0.,   0.,   0., 255., 255., 255.,   0.,   0.,   0.,
          255., 255.,],
 [   0., 255., 255., 255.,   0.,   0.,   0., 255., 255., 255.,   0.,   0.,   0.,
          255., 255.,],
 [ 255.,   0.,   0.,   0., 255., 255., 255.,   0.,   0.,   0., 255., 255., 255.,
            0.,   0.,],
 [ 255.,   0.,   0.,   0., 255., 255., 255.,   0.,   0.,   0., 255., 255., 255.,
            0.,   0.,],
 [ 255.,   0.,   0.,   0., 255., 255., 255.,   0.,   0.,   0., 255., 255., 255.,
            0.,   0.,]])


        image = FileImage(filename = resolve_path('topo/tests/unit/testimage.pgm'),
                      xdensity=15,
                      ydensity=15,
                      output_fns=[],
                      bounds=BoundingBox(radius=0.5))

        ps = image.pattern_sampler
        ps.size_normalization='fit_shortest'
        ps.whole_pattern_output_fns=[]

        assert_array_almost_equal(image_array,image())



    # CB: test rotation for PatternGenerators.
##     def test_rotation(self):
##         image_array = array(
## [[  96.59090909,  96.59090909,  96.59090909,  96.59090909,  96.59090909,
##          96.59090909,  96.59090909,  96.59090909,],
##  [   0.        ,  34.        ,  68.        , 102.        , 136.        ,
##         255.   ,   0.        ,   0.             ,],
##  [   0.        ,  34.        ,  68.        , 102.        , 136.        ,
##         255.        , 255.        ,   0.        ,],
##  [   0.        ,  34.        ,  68.        , 102.        , 136.        ,
##         255.        , 255.        , 255.        ,],
##  [ 255.        ,   0.        , 255.        ,   0.        , 255.        ,
##           0.        , 255.        ,   0.        ,],
##  [   0.        , 255.        ,   0.        , 255.        ,   0.        ,
##         255.        ,   0.        , 255.        ,],
##  [  96.59090909,  96.59090909,  96.59090909,  96.59090909,  96.59090909,
##          96.59090909,  96.59090909,  96.59090909,],
##  [  96.59090909,  96.59090909,  96.59090909,  96.59090909,  96.59090909,
##          96.59090909,  96.59090909,  96.59090909,]],Float)

##         image = FileImage(filename = resolve_path('topo/tests/unit/testimage.pgm'),
##                       density=8,
##                       bounds=BoundingBox(radius=0.5),
##                       output_fn=IdentityTF(),
##                       whole_image_output_fn=IdentityTF(),
##                       size_normalization='original')

##         assert_array_almost_equal(rot90(image_array,1),image(orientation=pi/2))
##         assert_array_almost_equal(rot90(image_array,-1),image(orientation=-pi/2))

if __name__ == "__main__":
	import nose
	nose.runmodule()