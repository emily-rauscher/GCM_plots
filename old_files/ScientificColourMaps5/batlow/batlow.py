# 
#         batlow
#                   www.fabiocrameri.ch/visualisation
from matplotlib.colors import LinearSegmentedColormap      
      
cm_data = [[0.0035963, 0.099991, 0.3501],      
           [0.0059481, 0.10412, 0.35081],      
           [0.0082407, 0.10821, 0.3515],      
           [0.010454, 0.11235, 0.3522],      
           [0.012856, 0.11642, 0.35291],      
           [0.014995, 0.12044, 0.3536],      
           [0.01708, 0.12449, 0.3543],      
           [0.019118, 0.12854, 0.355],      
           [0.021113, 0.13257, 0.35569],      
           [0.023067, 0.13656, 0.35638],      
           [0.024983, 0.14061, 0.35708],      
           [0.026865, 0.14462, 0.35777],      
           [0.028714, 0.14859, 0.35844],      
           [0.030534, 0.1526, 0.35913],      
           [0.032329, 0.15664, 0.35983],      
           [0.034083, 0.16061, 0.36051],      
           [0.036051, 0.16464, 0.3612],      
           [0.037779, 0.16866, 0.3619],      
           [0.039498, 0.17266, 0.36258],      
           [0.041209, 0.17667, 0.36326],      
           [0.042752, 0.1807, 0.36394],      
           [0.044421, 0.18476, 0.36463],      
           [0.045902, 0.1888, 0.36531],      
           [0.047477, 0.19284, 0.366],      
           [0.049077, 0.19691, 0.36668],      
           [0.050475, 0.20094, 0.36735],      
           [0.051912, 0.20504, 0.36803],      
           [0.053357, 0.20912, 0.36871],      
           [0.054819, 0.2132, 0.36938],      
           [0.056316, 0.21733, 0.37004],      
           [0.057649, 0.22144, 0.3707],      
           [0.059196, 0.22557, 0.37135],      
           [0.060593, 0.22968, 0.372],      
           [0.062035, 0.23383, 0.37264],      
           [0.063527, 0.23799, 0.37327],      
           [0.064987, 0.24215, 0.3739],      
           [0.066494, 0.24631, 0.37451],      
           [0.068052, 0.25049, 0.37509],      
           [0.069648, 0.2547, 0.37567],      
           [0.071235, 0.25888, 0.37624],      
           [0.072855, 0.26307, 0.37679],      
           [0.07453, 0.26725, 0.3773],      
           [0.076235, 0.27146, 0.37779],      
           [0.078026, 0.27569, 0.37825],      
           [0.079917, 0.27988, 0.37869],      
           [0.081925, 0.28408, 0.37909],      
           [0.083913, 0.28828, 0.37945],      
           [0.08592, 0.29248, 0.37978],      
           [0.088158, 0.29668, 0.38006],      
           [0.09042, 0.30085, 0.38029],      
           [0.092737, 0.30503, 0.38048],      
           [0.095262, 0.3092, 0.3806],      
           [0.097825, 0.31333, 0.38067],      
           [0.10046, 0.31745, 0.38068],      
           [0.10329, 0.32154, 0.38061],      
           [0.10621, 0.32562, 0.38048],      
           [0.10923, 0.32969, 0.38027],      
           [0.11238, 0.33371, 0.37998],      
           [0.11564, 0.33769, 0.37961],      
           [0.11905, 0.34165, 0.37916],      
           [0.12253, 0.34556, 0.37861],      
           [0.12623, 0.34945, 0.37798],      
           [0.13003, 0.35328, 0.37725],      
           [0.13391, 0.35707, 0.37643],      
           [0.13795, 0.36079, 0.37549],      
           [0.14208, 0.36449, 0.37447],      
           [0.14633, 0.36813, 0.37332],      
           [0.15071, 0.3717, 0.37208],      
           [0.1552, 0.37522, 0.37074],      
           [0.15977, 0.37868, 0.3693],      
           [0.1645, 0.38208, 0.36774],      
           [0.1693, 0.38542, 0.36609],      
           [0.17416, 0.3887, 0.36432],      
           [0.17913, 0.39191, 0.36246],      
           [0.18418, 0.39505, 0.36048],      
           [0.18931, 0.39814, 0.35843],      
           [0.19453, 0.40116, 0.35628],      
           [0.19975, 0.4041, 0.35402],      
           [0.20509, 0.40699, 0.35168],      
           [0.21047, 0.40981, 0.34927],      
           [0.21592, 0.41258, 0.34676],      
           [0.22139, 0.41527, 0.34417],      
           [0.22691, 0.41791, 0.34152],      
           [0.23246, 0.42048, 0.3388],      
           [0.23805, 0.42301, 0.336],      
           [0.24366, 0.42549, 0.33316],      
           [0.24933, 0.42791, 0.33025],      
           [0.25503, 0.43027, 0.32728],      
           [0.26071, 0.4326, 0.32427],      
           [0.26644, 0.43487, 0.32123],      
           [0.27216, 0.4371, 0.31814],      
           [0.27794, 0.43929, 0.31501],      
           [0.2837, 0.44145, 0.31185],      
           [0.2895, 0.44356, 0.30869],      
           [0.29528, 0.44564, 0.30549],      
           [0.30109, 0.44769, 0.30224],      
           [0.30693, 0.44972, 0.29901],      
           [0.31275, 0.45172, 0.29574],      
           [0.3186, 0.45369, 0.29247],      
           [0.32445, 0.45564, 0.2892],      
           [0.33033, 0.45757, 0.28589],      
           [0.3362, 0.45948, 0.28262],      
           [0.34209, 0.46138, 0.27932],      
           [0.34799, 0.46324, 0.27603],      
           [0.35389, 0.46511, 0.27271],      
           [0.35983, 0.46697, 0.26943],      
           [0.36577, 0.46881, 0.26614],      
           [0.37171, 0.47064, 0.26284],      
           [0.37769, 0.47246, 0.25957],      
           [0.38368, 0.47428, 0.25628],      
           [0.38968, 0.47609, 0.25304],      
           [0.39573, 0.47789, 0.24977],      
           [0.40177, 0.4797, 0.24655],      
           [0.40785, 0.48149, 0.24332],      
           [0.41396, 0.48329, 0.24012],      
           [0.42009, 0.48508, 0.23697],      
           [0.42626, 0.48688, 0.23378],      
           [0.43247, 0.48868, 0.23066],      
           [0.4387, 0.49048, 0.22756],      
           [0.44499, 0.49227, 0.22446],      
           [0.45132, 0.49408, 0.22145],      
           [0.45768, 0.49588, 0.21846],      
           [0.46412, 0.4977, 0.2155],      
           [0.47059, 0.4995, 0.2126],      
           [0.47712, 0.50134, 0.20975],      
           [0.48372, 0.50315, 0.20699],      
           [0.49039, 0.50499, 0.20428],      
           [0.49712, 0.50682, 0.20161],      
           [0.50393, 0.50867, 0.19906],      
           [0.5108, 0.51051, 0.19664],      
           [0.51777, 0.51237, 0.1943],      
           [0.5248, 0.51424, 0.19204],      
           [0.53193, 0.51609, 0.18994],      
           [0.53913, 0.51798, 0.18801],      
           [0.54644, 0.51984, 0.18619],      
           [0.55381, 0.52172, 0.18453],      
           [0.56129, 0.52361, 0.18305],      
           [0.56885, 0.52549, 0.18177],      
           [0.57651, 0.52738, 0.18071],      
           [0.58424, 0.52927, 0.17987],      
           [0.59207, 0.53114, 0.17926],      
           [0.59997, 0.53302, 0.1789],      
           [0.60795, 0.53488, 0.17878],      
           [0.616, 0.53675, 0.17893],      
           [0.62413, 0.5386, 0.17934],      
           [0.6323, 0.54044, 0.18004],      
           [0.64054, 0.54227, 0.18102],      
           [0.64882, 0.54408, 0.1823],      
           [0.65715, 0.54588, 0.18389],      
           [0.6655, 0.54765, 0.1858],      
           [0.67388, 0.54942, 0.18797],      
           [0.68228, 0.55115, 0.19039],      
           [0.69068, 0.55287, 0.19314],      
           [0.69909, 0.55457, 0.19617],      
           [0.70748, 0.55624, 0.19944],      
           [0.71587, 0.5579, 0.20302],      
           [0.72422, 0.55953, 0.20686],      
           [0.73255, 0.56114, 0.2109],      
           [0.74083, 0.56275, 0.21519],      
           [0.74907, 0.56432, 0.21974],      
           [0.75726, 0.56587, 0.22446],      
           [0.7654, 0.56742, 0.22943],      
           [0.77347, 0.56895, 0.23462],      
           [0.78147, 0.57047, 0.23995],      
           [0.7894, 0.57199, 0.24548],      
           [0.79725, 0.57351, 0.25122],      
           [0.80502, 0.57502, 0.25709],      
           [0.8127, 0.57654, 0.26314],      
           [0.82029, 0.57805, 0.26935],      
           [0.82778, 0.57959, 0.27571],      
           [0.83517, 0.58113, 0.28221],      
           [0.84245, 0.58268, 0.28885],      
           [0.84962, 0.58427, 0.29563],      
           [0.85666, 0.58588, 0.30254],      
           [0.86358, 0.5875, 0.30962],      
           [0.87037, 0.58915, 0.31679],      
           [0.87702, 0.59085, 0.32408],      
           [0.88353, 0.59259, 0.3315],      
           [0.88988, 0.59435, 0.33905],      
           [0.89607, 0.59615, 0.34668],      
           [0.90209, 0.598, 0.35442],      
           [0.90794, 0.59989, 0.36226],      
           [0.91361, 0.60183, 0.37019],      
           [0.91909, 0.60382, 0.37821],      
           [0.92437, 0.60584, 0.38629],      
           [0.92946, 0.60792, 0.39446],      
           [0.93433, 0.61004, 0.40269],      
           [0.93901, 0.6122, 0.41098],      
           [0.94346, 0.61441, 0.41931],      
           [0.9477, 0.61668, 0.42769],      
           [0.95172, 0.61898, 0.43609],      
           [0.95553, 0.62131, 0.44453],      
           [0.95911, 0.62369, 0.45297],      
           [0.96247, 0.6261, 0.46144],      
           [0.96562, 0.62855, 0.46991],      
           [0.96856, 0.63103, 0.47837],      
           [0.97129, 0.63354, 0.48684],      
           [0.97382, 0.63607, 0.4953],      
           [0.97616, 0.63863, 0.50374],      
           [0.9783, 0.64121, 0.51217],      
           [0.98026, 0.64382, 0.52058],      
           [0.98205, 0.64645, 0.52898],      
           [0.98367, 0.64909, 0.53735],      
           [0.98513, 0.65175, 0.54569],      
           [0.98644, 0.65441, 0.55401],      
           [0.98761, 0.6571, 0.56232],      
           [0.98865, 0.6598, 0.57059],      
           [0.98957, 0.66251, 0.57886],      
           [0.99037, 0.66521, 0.58709],      
           [0.99106, 0.66794, 0.59532],      
           [0.99166, 0.67068, 0.60353],      
           [0.99216, 0.67341, 0.61172],      
           [0.99259, 0.67616, 0.61991],      
           [0.99294, 0.67891, 0.62808],      
           [0.99321, 0.68166, 0.63624],      
           [0.99343, 0.68443, 0.6444],      
           [0.99359, 0.6872, 0.65256],      
           [0.9937, 0.68998, 0.66071],      
           [0.99376, 0.69277, 0.66888],      
           [0.99378, 0.69555, 0.67704],      
           [0.99376, 0.69834, 0.6852],      
           [0.99371, 0.70115, 0.69338],      
           [0.99363, 0.70395, 0.70156],      
           [0.99353, 0.70677, 0.70975],      
           [0.9934, 0.70959, 0.71796],      
           [0.99324, 0.71241, 0.72617],      
           [0.99307, 0.71524, 0.7344],      
           [0.99288, 0.71808, 0.74265],      
           [0.99267, 0.72092, 0.75091],      
           [0.99245, 0.72377, 0.75918],      
           [0.99222, 0.72663, 0.76747],      
           [0.99197, 0.7295, 0.77578],      
           [0.99171, 0.73237, 0.78411],      
           [0.99143, 0.73525, 0.79245],      
           [0.99115, 0.73813, 0.80082],      
           [0.99085, 0.74102, 0.8092],      
           [0.99055, 0.74392, 0.8176],      
           [0.99023, 0.74683, 0.82602],      
           [0.9899, 0.74973, 0.83446],      
           [0.98956, 0.75265, 0.84291],      
           [0.9892, 0.75558, 0.85138],      
           [0.98883, 0.75851, 0.85988],      
           [0.98844, 0.76144, 0.86839],      
           [0.98804, 0.76439, 0.87691],      
           [0.98762, 0.76733, 0.88546],      
           [0.98718, 0.77029, 0.89402],      
           [0.98673, 0.77325, 0.9026],      
           [0.98625, 0.77622, 0.9112],      
           [0.98575, 0.77919, 0.91981],      
           [0.98523, 0.78217, 0.92844],      
           [0.98468, 0.78515, 0.93708],      
           [0.9841, 0.78814, 0.94574],      
           [0.9835, 0.79113, 0.95441],      
           [0.98287, 0.79413, 0.9631],      
           [0.98221, 0.79713, 0.9718],      
           [0.98151, 0.80015, 0.98051]]      
      
batlow_map = LinearSegmentedColormap.from_list('batlow', cm_data)      
# For use of "viscm view"      
test_cm = batlow_map      
      
if __name__ == "__main__":      
    import matplotlib.pyplot as plt      
    import numpy as np      
      
    try:      
        from viscm import viscm      
        viscm(batlow_map)      
    except ImportError:      
        print("viscm not found, falling back on simple display")      
        plt.imshow(np.linspace(0, 100, 256)[None, :], aspect='auto',      
                   cmap=batlow_map)      
    plt.show()      
