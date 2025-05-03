import math
import matplotlib.pyplot as plt
from matplotlib.ticker import ScalarFormatter  

# from lp_large_n.txt (old version)
data = [
    [100000, 33642, 33646, 33668],
    [200000, 67703, 67704, 67739],
    [300000, 101903, 101906, 101945],
    [400000, 136143, 136149, 136186],
    [500000, 170456, 170459, 170502],
    [600000, 204811, 204815, 204858],
    [700000, 239187, 239196, 239241],
    [800000, 273604, 273610, 273668],
    [900000, 308029, 308042, 308099],
    [1000000, 342505, 342508, 342567],
    [2000000, 687796, 687800, 687883],
    [3000000, 1033949, 1033960, 1034056],
    [4000000, 1380625, 1380637, 1380747],
    [5000000, 1727605, 1727609, 1727731],
    [6000000, 2074962, 2074983, 2075114],
    [7000000, 2422486, 2422508, 2422651],
    [8000000, 2770212, 2770241, 2770389],
    [9000000, 3118129, 3118153, 3118302],
    [10000000, 3466235, 3466247, 3466414],
    [20000000, 6952243, 6952261, 6952477],
    [30000000, 10444441, 10444454, 10444694],
    [40000000, 13940484, 13940548, 13940838],
    [50000000, 17439282, 17439315, 17439638],
    [60000000, 20940210, 20940233, 20940591],
    [70000000, 24442818, 24442855, 24443233],
    [80000000, 27946958, 27947001, 27947403],
    [90000000, 31452431, 31452454, 31452859],
    [100000000, 34958725, 34959067, 34959207],
    [200000000, 70064782, 70065125, 70065426],
    [300000000, 105218403, 105218496, 105219155],
    [400000000, 140401212, 140401354, 140402099],
    [500000000, 175605266, 175605553, 175606238],
    [600000000, 210825848, 210825988, 210826906],
    [700000000, 246059851, 246060084, 246060998],
    [800000000, 281305291, 281306731, 281306449],
    [900000000, 316560601, 316560714, 316561839]
]

'''
t(10000) >= 3190 (heuristic fast) with (t-ceil(N/3)) = -144 (0.000s)
t(20000) >= 6525 (heuristic fast) with (t-ceil(N/3)) = -142 (0.000s)
t(30000) >= 9822 (heuristic fast) with (t-ceil(N/3)) = -178 (0.000s)
t(40000) >= 13169 (heuristic fast) with (t-ceil(N/3)) = -165 (0.001s)
t(50000) >= 16392 (heuristic fast) with (t-ceil(N/3)) = -275 (0.001s)
t(60000) >= 19690 (heuristic fast) with (t-ceil(N/3)) = -310 (0.001s)
t(70000) >= 23030 (heuristic fast) with (t-ceil(N/3)) = -304 (0.002s)
t(80000) >= 26530 (heuristic fast) with (t-ceil(N/3)) = -137 (0.001s)
t(90000) >= 29709 (heuristic fast) with (t-ceil(N/3)) = -291 (0.002s)
t(100000) >= 33184 (heuristic fast) with (t-ceil(N/3)) = -150 (0.001s)
t(200000) >= 66657 (heuristic fast) with (t-ceil(N/3)) = -10 (0.000s)
t(300000) >= 100408 (heuristic fast) with (t-ceil(N/3)) = 408 (0.000s)
t(400000) >= 133533 (heuristic fast) with (t-ceil(N/3)) = 199 (0.001s)
t(500000) >= 167365 (heuristic fast) with (t-ceil(N/3)) = 698 (0.000s)
t(600000) >= 201933 (heuristic fast) with (t-ceil(N/3)) = 1933 (0.000s)
t(700000) >= 235353 (heuristic fast) with (t-ceil(N/3)) = 2019 (0.001s)
t(800000) >= 269388 (heuristic fast) with (t-ceil(N/3)) = 2721 (0.001s)
t(900000) >= 303585 (heuristic fast) with (t-ceil(N/3)) = 3585 (0.001s)
t(1000000) >= 337642 (heuristic fast) with (t-ceil(N/3)) = 4308 (0.001s)
t(2000000) >= 679946 (heuristic fast) with (t-ceil(N/3)) = 13279 (0.001s)
t(3000000) >= 1023554 (heuristic fast) with (t-ceil(N/3)) = 23554 (0.002s)
t(4000000) >= 1365988 (heuristic fast) with (t-ceil(N/3)) = 32654 (0.003s)
t(5000000) >= 1706768 (heuristic fast) with (t-ceil(N/3)) = 40101 (0.003s)
t(6000000) >= 2049090 (heuristic fast) with (t-ceil(N/3)) = 49090 (0.003s)
t(7000000) >= 2390973 (heuristic fast) with (t-ceil(N/3)) = 57639 (0.004s)
t(8000000) >= 2738811 (heuristic fast) with (t-ceil(N/3)) = 72144 (0.005s)
t(9000000) >= 3082661 (heuristic fast) with (t-ceil(N/3)) = 82661 (0.004s)
t(10000000) >= 3422706 (heuristic fast) with (t-ceil(N/3)) = 89372 (0.004s)
t(20000000) >= 6879140 (heuristic fast) with (t-ceil(N/3)) = 212473 (0.008s)
t(30000000) >= 10334304 (heuristic fast) with (t-ceil(N/3)) = 334304 (0.007s)
t(40000000) >= 13833378 (heuristic fast) with (t-ceil(N/3)) = 500044 (0.012s)
t(50000000) >= 17252964 (heuristic fast) with (t-ceil(N/3)) = 586297 (0.013s)
t(60000000) >= 20676422 (heuristic fast) with (t-ceil(N/3)) = 676422 (0.021s)
t(70000000) >= 24156475 (heuristic fast) with (t-ceil(N/3)) = 823141 (0.019s)
t(80000000) >= 27691895 (heuristic fast) with (t-ceil(N/3)) = 1025228 (0.021s)
t(90000000) >= 31116792 (heuristic fast) with (t-ceil(N/3)) = 1116792 (0.025s)
t(100000000) >= 34616026 (heuristic fast) with (t-ceil(N/3)) = 1282692 (0.021s)
t(200000000) >= 69326602 (heuristic fast) with (t-ceil(N/3)) = 2659935 (0.044s)
t(300000000) >= 104262400 (heuristic fast) with (t-ceil(N/3)) = 4262400 (0.048s)
t(400000000) >= 139136498 (heuristic fast) with (t-ceil(N/3)) = 5803164 (0.052s)
t(500000000) >= 173960145 (heuristic fast) with (t-ceil(N/3)) = 7293478 (0.081s)
t(600000000) >= 209024651 (heuristic fast) with (t-ceil(N/3)) = 9024651 (0.111s)
t(700000000) >= 244134457 (heuristic fast) with (t-ceil(N/3)) = 10801123 (0.094s)
t(800000000) >= 278817959 (heuristic fast) with (t-ceil(N/3)) = 12151292 (0.108s)
t(900000000) >= 313422748 (heuristic fast) with (t-ceil(N/3)) = 13422748 (0.111s)
t(1000000000) >= 348659617 (heuristic fast) with (t-ceil(N/3)) = 15326283 (0.132s)
t(2000000000) >= 699375340 (heuristic fast) with (t-ceil(N/3)) = 32708673 (0.195s)
t(3000000000) >= 1049265763 (heuristic fast) with (t-ceil(N/3)) = 49265763 (0.246s)
t(4000000000) >= 1399130424 (heuristic fast) with (t-ceil(N/3)) = 65797090 (0.384s)
t(5000000000) >= 1752451206 (heuristic fast) with (t-ceil(N/3)) = 85784539 (0.474s)
t(6000000000) >= 2102462426 (heuristic fast) with (t-ceil(N/3)) = 102462426 (0.432s)
t(7000000000) >= 2452155258 (heuristic fast) with (t-ceil(N/3)) = 118821924 (0.476s)
t(8000000000) >= 2802515001 (heuristic fast) with (t-ceil(N/3)) = 135848334 (0.794s)
t(9000000000) >= 3156503781 (heuristic fast) with (t-ceil(N/3)) = 156503781 (0.584s)
t(10000000000) >= 3506181640 (heuristic fast) with (t-ceil(N/3)) = 172848306 (0.593s)
t(20000000000) >= 7023289176 (heuristic fast) with (t-ceil(N/3)) = 356622509 (1.382s)
t(30000000000) >= 10550213184 (heuristic fast) with (t-ceil(N/3)) = 550213184 (1.559s)
t(40000000000) >= 14057249297 (heuristic fast) with (t-ceil(N/3)) = 723915963 (2.255s)
t(50000000000) >= 17591426970 (heuristic fast) with (t-ceil(N/3)) = 924760303 (2.745s)
t(60000000000) >= 21114347590 (heuristic fast) with (t-ceil(N/3)) = 1114347590 (2.289s)
t(70000000000) >= 24654219219 (heuristic fast) with (t-ceil(N/3)) = 1320885885 (2.954s)
t(80000000000) >= 28166354328 (heuristic fast) with (t-ceil(N/3)) = 1499687661 (2.303s)
t(90000000000) >= 31699844022 (heuristic fast) with (t-ceil(N/3)) = 1699844022 (3.292s)
t(100000000000) >= 35207945082 (heuristic fast) with (t-ceil(N/3)) = 1874611748 (3.028s)
t(200000000000) >= 70582291950 (heuristic fast) with (t-ceil(N/3)) = 3915625283 (5.544s)
t(300000000000) >= 105839929827 (heuristic fast) with (t-ceil(N/3)) = 5839929827 (8.476s)
t(400000000000) >= 141202218034 (heuristic fast) with (t-ceil(N/3)) = 7868884700 (10.238s)
t(500000000000) >= 176660075193 (heuristic fast) with (t-ceil(N/3)) = 9993408526 (8.887s)
t(600000000000) >= 212000894709 (heuristic fast) with (t-ceil(N/3)) = 12000894709 (13.963s)
t(700000000000) >= 247600693913 (heuristic fast) with (t-ceil(N/3)) = 14267360579 (15.827s)
t(800000000000) >= 282759334392 (heuristic fast) with (t-ceil(N/3)) = 16092667725 (13.326s)
t(900000000000) >= 318000733240 (heuristic fast) with (t-ceil(N/3)) = 18000733240 (18.559s)
t(1000000000000) >= 353566331002 (heuristic fast) with (t-ceil(N/3)) = 20232997668 (16.917s)'''
# from tbounds.txt, renamed to tbounds_heuristic_fast_greedy_1e4_1e12.txt
data2 = [
    [10000, 3190, 0.000],
    [20000, 6525, 0.000],
    [30000, 9822, 0.000],
    [40000, 13169, 0.001],
    [50000, 16392, 0.001],
    [60000, 19690, 0.001],
    [70000, 23030, 0.002],
    [80000, 26530, 0.001],
    [90000, 29709, 0.002],
    [100000, 33184, 0.001],
    [200000, 66657, 0.000],
    [300000, 100408, 0.000],
    [400000, 133533, 0.001],
    [500000, 167365, 0.000],
    [600000, 201933, 0.000],
    [700000, 235353, 0.001],   
    [800000, 269388, 0.001],
    [900000, 303585, 0.001],
    [1000000, 337642, 0.001],
    [2000000, 679946, 0.001],
    [3000000, 1023554, 0.002],
    [4000000, 1365988, 0.003],
    [5000000, 1706768, 0.003],
    [6000000, 2049090, 0.003],
    [7000000, 2390973, 0.004],
    [8000000, 2738811, 0.005],
    [9000000, 3082661, 0.004],
    [10000000, 3422706, 0.004],
    [20000000, 6879140, 0.008],
    [30000000, 10334304, 0.007],
    [40000000, 13833378, 0.012],
    [50000000, 17252964, 0.013],
    [60000000, 20676422, 0.021],
    [70000000, 24156475, 0.019],
    [80000000, 27691895, 0.021],
    [90000000, 31116792, 0.025],
    [100000000, 34616026, 0.021],
    [200000000, 69326602, 0.044],
    [300000000, 104262400, 0.048],
    [400000000, 139136498, 0.052],
    [500000000, 173960145, 0.081],
    [600000000, 209024651, 0.111],
    [700000000, 244134457, 0.094],
    [800000000, 278817959, 0.108],
    [900000000, 313422748, 0.111],
    [1000000000, 348659617, 0.132],
    [2000000000, 699375340, 0.195],
    [3000000000, 1049265763, 0.246],
    [4000000000, 1399130424, 0.384],
    [5000000000, 1752451206, 0.474],
    [6000000000, 2102462426, 0.432],
    [7000000000, 2452155258, 0.476],
    [8000000000, 2802515001, 0.794],
    [9000000000, 3156503781, 0.584],
    [10000000000, 3506181640, 0.593],
    [20000000000, 7023289176, 1.382],
    [30000000000, 10550213184, 1.559],
    [40000000000, 14057249297, 2.255],
    [50000000000, 17591426970, 2.745],
    [60000000000, 21114347590, 2.289],
    [70000000000, 24654219219, 2.954],
    [80000000000, 28166354328, 2.303],
    [90000000000, 31699844022, 3.292],
    [100000000000,35207945082 ,3.028],
    [200000000000, 70582291950, 5.544],
    [300000000000, 105839929827, 8.476],
    [400000000000, 141202218034, 10.238],
    [500000000000, 176660075193, 8.887],
    [600000000000, 212000894709, 13.963],
    [700000000000, 247600693913, 15.827],
    [800000000000, 282759334392, 13.326],
    [900000000000, 318000733240, 18.559],
    [1000000000000, 353566331002, 16.917]    
]

# from tbounds_o.txt, renamed to tbounds_exhaustive_fast_greedy_1e4_1e8.txt
data3 = [
    [10000, 3190, 0.008],
    [20000, 6525, 0.003],
    [30000, 9822, 0.001],
    [40000, 13169, 0.001],
    [50000, 16392, 0.002],
    [60000, 19690, 0.002],
    [70000, 23030, 0.002],
    [80000, 26530, 0.002],
    [90000, 29709, 0.003],
    [100000, 33184, 0.003],
    [200000, 66657, 0.003],
    [300000, 100408, 0.006],
    [400000, 134536, 0.011],
    [500000, 168213, 0.016],
    [600000, 202612, 0.021],
    [700000, 236398, 0.037],
    [800000, 270220, 0.033],
    [900000, 304032, 0.040],
    [1000000, 338782, 0.050],
    [2000000, 680922, 0.112],
    [3000000, 1024158, 0.180],
    [4000000, 1369123, 0.294],
    [5000000, 1711827, 0.515],
    [6000000, 2053781, 0.668],
    [7000000, 2400144, 0.882],
    [8000000, 2742579, 0.988],
    [9000000, 3088060, 1.179],
    [10000000, 3433084, 1.586],
    [20000000, 6889761, 4.258],
    [30000000, 10346651, 8.120],
    [40000000, 13835949, 9.294],
    [50000000, 17291487, 19.873],
    [60000000, 20711660, 31.423],
    [70000000, 24208875, 37.419],
    [80000000, 27733589, 35.919],
    [90000000, 31176678, 51.054],
    [100000000, 34673154, 55.608]
]


# from lp_large_n.txt
data4 = [
    [ 10000, 3252, 3269, 3271, 3278],
    [ 20000, 6588, 6604, 6605, 6616],
    [  30000, 9943, 9962, 9964, 9977],
    [ 40000, 13304, 13325, 13327, 13342],
    [ 50000, 16669, 16700, 16701, 16716],
    [ 60000, 20049, 20074, 20077, 20094],
    [ 70000, 23433, 23459, 23461, 23483],
    [ 80000, 26829, 26856, 26858, 26875],
    [ 90000, 30213, 30239, 30244, 30266],
    [ 100000, 33614, 33642, 33646, 33668],
    [ 200000, 67667, 67703, 67704, 67739],
    [ 300000, 101854, 101903, 101906, 101945],
    [ 400000, 136090, 136143, 136149, 136186],
    [ 500000, 170402, 170456, 170459, 170502],
    [ 600000, 204752, 204811, 204815, 204858],
    [ 700000, 239135, 239187, 239196, 239241],
    [ 800000, 273568, 273604, 273610, 273668],
    [ 900000, 307968, 308029, 308042, 308099],
    [1000000, 342450, 342505, 342508, 342567],
    [2000000, 687735, 687796, 687800, 687883],
    [3000000, 1033857, 1033949, 1033960, 1034056],
    [4000000, 1380531, 1380625, 1380637, 1380747],
    [5000000, 1727499, 1727605, 1727609, 1727731],
    [6000000, 2074830, 2074962, 2074983, 2075114],
    [7000000, 2422422, 2422486, 2422508, 2422651],
    [8000000, 2770115, 2770212, 2770241, 2770389],
    [9000000, 3118039, 3118129, 3118153, 3118302],
    [10000000, 3466108, 3466235, 3466247, 3466414],
    [20000000, 6952115, 6952243, 6952261, 6952477],
    [30000000, 10444300, 10444441, 10444454, 10444694],
    [40000000, 13940348, 13940484, 13940548, 13940838],
    [50000000, 17439056, 17439282, 17439315, 17439638],
    [60000000, 20940044, 20940210, 20940233, 20940591],
    [70000000, 24442600, 24442818, 24442855, 24443233],
    [80000000, 27946752, 27946958, 27947001, 27947403],
    [90000000, 31452204, 31452431, 31452454, 31452859],
    [100000000, 34958488, 34958725, 34958773, 34959207],
    [200000000, 70064542, 70064782, 70064827, 70065426],
    [300000000, 105218127, 105218403, 105218444, 105219155],
    [400000000, 140400896, 140401212, 140401292, 140402099],
    [500000000, 175604950, 175605266, 175605364, 175606238],
    [600000000, 210825475, 210825848, 210825916, 210826906],
    [700000000, 246059457, 246059851, 246059940, 246060998],
    [800000000, 281304919, 281305291, 281305383, 281306449],
    [900000000, 316560258, 316560601, 316560702, 316561839]
]

# from tbounds_g.txt, renamed to tbounds_exhaustive_greedy_1e4_1e8.txt
data5 = [
    [10000, 3258, 0.003],
    [20000, 6578, 0.002],
    [30000, 9912, 0.004],
    [40000, 13303, 0.006],
    [50000, 16667, 0.008],
    [60000, 19950, 0.005],
    [70000, 23414, 0.003],
    [80000, 26791, 0.007],
    [90000, 30177, 0.004],
    [100000, 33572, 0.015],
    [200000, 67613, 0.016],
    [300000, 101849, 0.020],
    [400000, 135996, 0.046],
    [500000, 170405, 0.067],
    [600000, 204597, 0.071],
    [700000, 238908, 0.066],
    [800000, 273378, 0.096],
    [900000, 307803, 0.190],
    [1000000, 342303, 0.257],
    [2000000, 687503, 0.370],
    [3000000, 1033385, 0.496],
    [4000000, 1379920, 0.731],
    [5000000, 1727135, 1.983],
    [6000000, 2073482, 2.745],
    [7000000, 2421657, 9.079],
    [8000000, 2769029, 5.865],
    [9000000, 3117029, 3.997],
    [10000000, 3465013, 4.091],
    [20000000, 6949513, 14.860],
    [30000000, 10442788, 55.389],
    [40000000, 13938940, 29.764],
    [50000000, 17437371, 130.826],
    [60000000, 20938423, 129.957],
    [70000000, 24441771, 528.612],
    [80000000, 27943125, 395.025],
    [90000000, 31448310, 1581.736],
    [100000000,34955940 ,613.914]
]


'''
t(10000) >= 3258 (heuristic greedy) with (t-ceil(N/3)) = -76 (0.000s)
t(20000) >= 6578 (heuristic greedy) with (t-ceil(N/3)) = -89 (0.002s)
t(30000) >= 9912 (heuristic greedy) with (t-ceil(N/3)) = -88 (0.003s)
t(40000) >= 13303 (heuristic greedy) with (t-ceil(N/3)) = -31 (0.001s)
t(50000) >= 16667 (heuristic greedy) with (t-ceil(N/3)) = 0 (0.000s)
t(60000) >= 19950 (heuristic greedy) with (t-ceil(N/3)) = -50 (0.004s)
t(70000) >= 23380 (heuristic greedy) with (t-ceil(N/3)) = 46 (0.001s)
t(80000) >= 26767 (heuristic greedy) with (t-ceil(N/3)) = 100 (0.001s)
t(90000) >= 30172 (heuristic greedy) with (t-ceil(N/3)) = 172 (0.001s)
t(100000) >= 33337 (heuristic greedy) with (t-ceil(N/3)) = 3 (0.002s)
t(200000) >= 67482 (heuristic greedy) with (t-ceil(N/3)) = 815 (0.003s)
t(300000) >= 101708 (heuristic greedy) with (t-ceil(N/3)) = 1708 (0.003s)
t(400000) >= 135728 (heuristic greedy) with (t-ceil(N/3)) = 2394 (0.005s)
t(500000) >= 169893 (heuristic greedy) with (t-ceil(N/3)) = 3226 (0.005s)
t(600000) >= 204191 (heuristic greedy) with (t-ceil(N/3)) = 4191 (0.006s)
t(700000) >= 238765 (heuristic greedy) with (t-ceil(N/3)) = 5431 (0.008s)
t(800000) >= 273035 (heuristic greedy) with (t-ceil(N/3)) = 6368 (0.010s)
t(900000) >= 306791 (heuristic greedy) with (t-ceil(N/3)) = 6791 (0.011s)
t(1000000) >= 340790 (heuristic greedy) with (t-ceil(N/3)) = 7456 (0.012s)
t(2000000) >= 686738 (heuristic greedy) with (t-ceil(N/3)) = 20071 (0.028s)
t(3000000) >= 1033006 (heuristic greedy) with (t-ceil(N/3)) = 33006 (0.042s)
t(4000000) >= 1379632 (heuristic greedy) with (t-ceil(N/3)) = 46298 (0.059s)
t(5000000) >= 1724823 (heuristic greedy) with (t-ceil(N/3)) = 58156 (0.092s)
t(6000000) >= 2071720 (heuristic greedy) with (t-ceil(N/3)) = 71720 (0.085s)
t(7000000) >= 2412529 (heuristic greedy) with (t-ceil(N/3)) = 79195 (0.118s)
t(8000000) >= 2764900 (heuristic greedy) with (t-ceil(N/3)) = 98233 (0.135s)
t(9000000) >= 3115006 (heuristic greedy) with (t-ceil(N/3)) = 115006 (0.176s)
t(10000000) >= 3463360 (heuristic greedy) with (t-ceil(N/3)) = 130026 (0.179s)
t(20000000) >= 6946752 (heuristic greedy) with (t-ceil(N/3)) = 280085 (0.444s)
t(30000000) >= 10429746 (heuristic greedy) with (t-ceil(N/3)) = 429746 (0.616s)
t(40000000) >= 13935056 (heuristic greedy) with (t-ceil(N/3)) = 601722 (0.772s)
t(50000000) >= 17421954 (heuristic greedy) with (t-ceil(N/3)) = 755287 (1.091s)
t(60000000) >= 20927883 (heuristic greedy) with (t-ceil(N/3)) = 927883 (1.089s)
t(70000000) >= 24401463 (heuristic greedy) with (t-ceil(N/3)) = 1068129 (1.638s)
t(80000000) >= 27918926 (heuristic greedy) with (t-ceil(N/3)) = 1252259 (1.551s)
t(90000000) >= 31354920 (heuristic greedy) with (t-ceil(N/3)) = 1354920 (1.538s)
t(100000000) >= 34924234 (heuristic greedy) with (t-ceil(N/3)) = 1590900 (1.826s)
t(200000000) >= 69835925 (heuristic greedy) with (t-ceil(N/3)) = 3169258 (3.354s)
t(300000000) >= 105133902 (heuristic greedy) with (t-ceil(N/3)) = 5133902 (4.439s)
t(400000000) >= 140280238 (heuristic greedy) with (t-ceil(N/3)) = 6946904 (7.564s)
t(500000000) >= 175358640 (heuristic greedy) with (t-ceil(N/3)) = 8691973 (10.269s)
t(600000000) >= 210712478 (heuristic greedy) with (t-ceil(N/3)) = 10712478 (11.349s)
t(700000000) >= 245882180 (heuristic greedy) with (t-ceil(N/3)) = 12548846 (14.494s)
t(800000000) >= 281199770 (heuristic greedy) with (t-ceil(N/3)) = 14533103 (17.530s)
t(900000000) >= 316308547 (heuristic greedy) with (t-ceil(N/3)) = 16308547 (22.831s)
t(1000000000) >= 351661732 (heuristic greedy) with (t-ceil(N/3)) = 18328398 (16.854s)'''
# from tbounds_heuristic_greedy_1e4_1e9.txt
data6 = [
    [10000, 3258, 0.000],
    [20000, 6578, 0.002],
    [30000, 9912, 0.003],
    [40000, 13303, 0.001],
    [50000, 16667, 0.000],
    [60000, 19950, 0.004],
    [70000, 23380, 0.001],
    [80000, 26767, 0.001],
    [90000, 30172, 0.001],
    [100000, 33337, 0.002],
    [200000, 67482, 0.003],
    [300000, 101708, 0.003],
    [400000, 135728, 0.005],
    [500000, 169893, 0.005],
    [600000, 204191, 0.006],
    [700000, 238765, 0.008],
    [800000, 273035, 0.010],
    [900000, 306791, 0.011],
    [1000000, 340790, 0.012],
    [2000000, 686738, 0.028],
    [3000000, 1033006, 0.042],
    [4000000, 1379632, 0.059],
    [5000000, 1724823, 0.092],
    [6000000, 2071720, 0.085],
    [7000000, 2412529, 0.118],
    [8000000, 2764900, 0.135],
    [9000000, 3115006, 0.176],
    [10000000,3463360 ,0.179],
    [20000000, 6946752, 0.444],
    [30000000, 10429746, 0.616],
    [40000000, 13935056, 0.772],
    [50000000, 17421954, 1.091],
    [60000000, 20927883, 1.089],
    [70000000, 24401463, 1.638],
    [80000000, 27918926, 1.551],
    [90000000, 31354920, 1.538],
    [100000000,34924234 ,1.826],
    [200000000, 69835925, 3.354],
    [300000000, 105133902, 4.439],
    [400000000, 140280238, 7.564],
    [500000000, 175358640, 10.269],
    [600000000, 210712478, 11.349],
    [700000000, 245882180, 14.494],
    [800000000, 281199770, 17.530],
    [900000000, 316308547, 22.831],
    [1000000000,351661732 ,16.854]
]

def is_prime(n):
    """Return True if n is a prime number, else False."""
    if n < 2:
        return False
    # Check divisibility from 2 up to sqrt(n)
    for i in range(2, int(math.sqrt(n)) + 1):
        if n % i == 0:
            return False
    return True

def logfac(N):
    sum = 0
    for i in range(1, N+1):
        sum += math.log(i)
    return sum

# Test if  \sum_{p > \frac{t}{\sqrt{t}+1}} \left\lfloor \frac{N}{p} \right\rfloor \log \left( \frac{p}{t} \left\lceil \frac{t}{p} \right\rceil \right) > \log N! - N \log t

def criterion( N, t ):
    sum = 0
    threshold = logfac(N) - N * math.log(t) + 0.0001 
    tfail = t/(math.floor(math.sqrt(t)))
    for p in reversed(range(1, N + 1)):
        if p < tfail:
            return False
        if is_prime(p):
            sum += math.floor(N/p) * math.log((p/t)*math.ceil(t/p))
            if sum > threshold:
                return True
    return False

def best_t( N, init ):
    t = init
    while not criterion(N, t):
        t += 10
    while criterion(N,t):
        t -= 1
    print(f"Testing {N}: needed to increment by {t-init}") 
    return t



c0 = 0.304419010
c1 = 0.75554808

def round_to_sigfigs(n, sigfigs):
    if n == 0:
        return 0
    import math
    digits = int(math.floor(math.log10(abs(n)))) + 1
    factor = 10**(max(digits - sigfigs,0))
    return (n // factor) * factor

def table1():
    for x in data:
        N, t_lower, t_upper, lemma = x
        diff = t_upper - t_lower
        lemmadiff = lemma - t_lower
        approx = int(N / math.e - c0 * N / math.log(N) - c1 * N / (math.log(N) ** 2))
        approxdiff = t_lower - approx  
        digits = int(math.floor(math.log10(N))) 

        print(f"${N//10**digits} \\times 10^{digits}$ & $\\num{{{t_lower}}}$ & $t(N)_- + {diff}$ & $t(N)_- + {lemmadiff}$ & $t(N)_- - \\num{{{approxdiff}}}$ \\\\")

def table2():
    for x in data:
        N, t_lower, _, _ = x
        t_greedy = t_lower
        time_greedy = 0
        t_greedy_2 = t_lower
        time_greedy_2 = 0

        for y in data2:
            if y[0] == N:
                t_greedy = y[1]
                time_greedy = y[2]
                break
        for y in data3:
            if y[0] == N:
                t_greedy_2 = y[1]
                time_greedy_2 = y[2]
                break

        digits = int(math.floor(math.log10(N))) 
        print(f"${N//10**digits} \\times 10^{digits}$ & $\\num{{{t_lower}}}$ & $t(N)_- - \\num{{{t_lower-t_greedy}}}$ & {time_greedy} & $t(N)_- - \\num{{{t_lower-t_greedy_2}}}$ & {time_greedy_2} \\\\")


def plot():
    base1 = []
    t1 = []
    for x in data:
        N, t_lower, _, _ = x
        base1.append(N)
        t1.append(t_lower / N)
    

    base2 = []
    t2 = []
    for x in data2:
        N, t, _ = x
        base2.append(N)
        t2.append(t / N)

    base3 = []
    t3 = []
    for x in data3:
        N, t, _ = x
        base3.append(N)
        t3.append(t / N)

    base5 = []
    t5 = []
    for x in data5:
        N, t, _ = x
        base5.append(N)
        t5.append(t / N)


    # form the union of bases
    base = sorted(set(base1) | set(base2) | set(base5))

    asym = [1/math.e - 0.30440119010/math.log(N) for N in base]
    asym2 = [1/math.e - 0.30440119010/math.log(N) - 0.75554808/math.log(N)**2 for N in base]

    fig, ax = plt.subplots()
    ax.plot(base, asym, linestyle="--", label='$1/e-c_0/\\log N$', color='purple' )
    ax.plot(base, asym2, linestyle="--", label='$1/e-c_0/\\log N-c_1/\\log^2 N$', color='gray' )
    ax.plot(base, [1/math.e for _ in base], linestyle="--", label='$1/e$', color='orange' )
    ax.plot(base, [1/3 for _ in base], linestyle="--", label='$1/3$', color='red' )
    ax.plot(base1, t1, label='LP lower bound', color='blue' )
    ax.plot(base2, t2, label='Heuristic fast greedy', color='green' )
    ax.plot(base3, t3, label='Exhaustive fast greedy', color='brown' )
    ax.plot(base5, t5, label='Exhaustive greedy', color='pink' )
    ax.set_xscale('log') 
#    ax.xaxis.set_major_formatter(ScalarFormatter())
#    ax.ticklabel_format(style='plain', axis='x')            # turn off offset
    plt.title('Approximations to $t(N)/N$')
    plt.xlabel('$N$ (log scale)')
    plt.ylim(0.33,0.37)
    plt.legend()
    plt.grid(True)
    plt.tight_layout()
    plt.show()


def plot2():
    base1 = []
    t1 = []
    for x in data:
        N, t_lower, _, _ = x
        base1.append(N)
        t1.append(t_lower / N)
    

    base2 = []
    t2 = []
    for x in data2:
        N, t, _ = x
        base2.append(N)
        t2.append(t / N)

    base3 = []
    t3 = []
    for x in data3:
        N, t, _ = x
        base3.append(N)
        t3.append(t / N)

    base4 = []
    t4_floor = []
    t4_lower = []
    t4_upper = []
    t4_bound = []
    for x in data4:
        N, t_floor, t_lower, t_upper, t_bound = x
        base4.append(N)
        t4_floor.append(t_floor / N)
        t4_lower.append(t_lower / N)
        t4_upper.append(t_upper / N)
        t4_bound.append(t_bound / N)

    base5 = []
    t_greedy_opt = []
    for x in data5:
        N, t, _ = x
        base5.append(N)
        t_greedy_opt.append(t / N)

    base6 = []
    t_greedy_heuristic = []
    for x in data6:
        N, t, _ = x
        base6.append(N)
        t_greedy_heuristic.append(t / N)


    # form the union of base1 and base2
    base = sorted(set(base1) | set(base2) | set(base3) | set(base4) | set(base5) | set(base6))
    asym = [1/math.e - 0.30440119010/math.log(N) for N in base]
    asym2 = [1/math.e - 0.30440119010/math.log(N) - 0.75554808/math.log(N)**2 for N in base]

    fig, ax = plt.subplots()
    ax.plot(base, asym, linestyle="--", label='$1/e-c_0/\\log N$', color='purple' )
    ax.plot(base, asym2, linestyle="--", label='$1/e-c_0/\\log N-c_1/\\log^2 N$', color='gray' )
#    ax.plot(base, [1/math.e for _ in base], linestyle="--", label='$1/e$', color='orange' )
    ax.plot(base, [1/3 for _ in base], linestyle="--", label='$1/3$', color='red' )
    ax.plot(base4, t4_lower, label='LP lower bound', color='blue' )
    ax.plot(base4, t4_upper, label='LP upper bound', color='cyan' )
    ax.plot(base2, t2, label='Heuristic fast greedy', color='green' )
    ax.plot(base3, t3, label='Exhaustive fast greedy', color='brown' )
    ax.plot(base4, t4_floor, label='LP floor bound', color='pink' )
    ax.plot(base4, t4_bound, label='Lemma 5.1', color='black' )
    ax.plot(base6, t_greedy_heuristic, label='Heuristic greedy', color='olive' )
    ax.plot(base5, t_greedy_opt, label='Exhaustive greedy', color='orange' )
    
    ax.set_xscale('log') 
#    ax.xaxis.set_major_formatter(ScalarFormatter())
#    ax.ticklabel_format(style='plain', axis='x')            # turn off offset
    plt.title('Approximations to $t(N)/N$')
    plt.xlabel('$N$ (log scale)')
    plt.ylim(0.33,0.355)
    plt.legend()
    plt.grid(True)
    plt.tight_layout()
    plt.show()

plot2()