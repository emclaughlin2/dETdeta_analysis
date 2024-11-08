#include <iostream>
#include <TH2D.h>
#include <TH1D.h>
#include <TChain.h>
#include <TMath.h>
#include <TTree.h>
#include <TFile.h>
#include <sstream> //std::ostringstsream
#include <fstream> //std::ifstream
#include <iostream> //std::cout, std::endl
#include <cmath>
#include <TGraphErrors.h>
#include <TGraph.h>
#include <TSpectrum.h>
#include <TF1.h>
#include <TTreeReader.h>
#include <TTreeReaderValue.h>
#include <TTreeReaderArray.h>
#include <string>
#include <set>
#include <TVector3.h>
#include <map>
#include <vector>
#include <TDatabasePDG.h>
#include <tuple>

using namespace std;

std::set<std::tuple<int, int>> emcal_hot_dead_map = {{77,186},{79,187}};
std::set<std::tuple<int, int>> ihcal_hot_dead_map = {{8,32},{7,51}};
std::set<std::tuple<int, int>> ohcal_hot_dead_map;

double emcal_eta_bin_centers[24];
double ihcal_eta_bin_centers[24];
double ohcal_eta_bin_centers[24];
double calo_eta_bin_centers[24];

double vertex_reweight[200] = {1.0};
std::vector<float> centrality_bin;

const double eta_bin_centers[24] = {-1.05417,-0.9625,-0.870833,-0.779167,-0.6875,-0.595833,-0.504167,-0.4125,-0.320833,-0.229167,
	-0.1375,-0.0458333,0.0458333,0.1375,0.229167,0.320833,0.4125,0.504167,0.595833,0.6875,0.779167,0.870833,0.9625,1.05417};

float em_zs_calib[96][256];
float ih_zs_calib[24][64];
float oh_zs_calib[24][64];

const int cemcSize = 24576;
const int ihcalSize = 1536;
const int ohcalSize = 1536;
const int mbdSize = 256;
const int g4Size = 100000;
const int vtxSize = 3;

std::set<long> skip_event = {3, 5, 81, 82, 83, 84, 86, 87, 88, 89, 90, 250, 302, 308, 316, 353, 391, 392, 398, 400, 446, 449, 452, 454, 455, 
532, 534, 536, 537, 538, 540, 583, 587, 588, 634, 635, 701, 787, 813, 815, 819, 824, 871, 875, 877, 880, 971, 973, 1046, 1061, 1099, 1131, 
1132, 1133, 1134, 1135, 1136, 1138, 1139, 1140, 1164, 1168, 1170, 1221, 1223, 1257, 1330, 1339, 1414, 1418, 1454, 1457, 1505, 1580, 1691, 1698, 
1721, 1730, 1881, 1882, 1886, 1887, 1902, 1942, 1944, 1945, 1948, 1949, 1950, 1953, 1954, 1955, 1956, 1960, 2001, 2008, 2009, 2036, 2062, 2069, 
2081, 2082, 2086, 2088, 2089, 2090, 2153, 2154, 2156, 2160, 2181, 2182, 2183, 2184, 2185, 2186, 2187, 2188, 2189, 2190, 2217, 2218, 2240, 2252, 
2253, 2255, 2256, 2258, 2260, 2271, 2273, 2274, 2277, 2278, 2279, 2302, 2304, 2306, 2307, 2308, 2314, 2319, 2398, 2431, 2432, 2433, 2434, 2493, 
2572, 2573, 2574, 2575, 2576, 2581, 2585, 2589, 2671, 2673, 2863, 2864, 2866, 2867, 2868, 2869, 2981, 2983, 2985, 2987, 2992, 2994, 2995, 2997, 
2998, 3013, 3015, 3018, 3019, 3367, 3435, 3553, 3557, 3565, 3566, 3567, 3569, 3602, 3605, 3606, 3609, 3660, 3671, 3673, 3674, 3677, 3680, 3692, 
3695, 3724, 3727, 3730, 3835, 3840, 3901, 3903, 3911, 3912, 3914, 3915, 3918, 3921, 3922, 3923, 3924, 3925, 3926, 3927, 3928, 3929, 3930, 3955, 
4034, 4046, 4048, 4062, 4063, 4065, 4068, 4069, 4091, 4093, 4094, 4099, 4221, 4222, 4224, 4225, 4228, 4229, 4230, 4241, 4244, 4248, 4361, 4362, 
4363, 4364, 4368, 4370, 4371, 4401, 4402, 4404, 4405, 4406, 4441, 4443, 4462, 4466, 4564, 4592, 4596, 4598, 4600, 4633, 4653, 4692, 4694, 4695, 
4698, 4734, 4735, 4736, 4737, 4738, 4763, 4771, 4773, 4776, 4777, 4811, 4812, 4813, 4816, 4818, 4862, 4863, 4864, 4865, 4866, 4868, 4901, 4902, 
4903, 4905, 4907, 4909, 4947, 4970, 4985, 4986, 4987, 4988, 4990, 5029, 5031, 5185, 5298, 5313, 5341, 5342, 5344, 5345, 5347, 5349, 5350, 5361, 
5362, 5366, 5367, 5368, 5371, 5373, 5377, 5380, 5469, 5516, 5517, 5520, 5583, 5591, 5592, 5593, 5594, 5595, 5596, 5597, 5621, 5763, 5764, 5769, 
5781, 5786, 5821, 5862, 5868, 5894, 5937, 5942, 5943, 5944, 5949, 6108, 6122, 6128, 6130, 6204, 6205, 6206, 6207, 6208, 6209, 6256, 6257, 6396, 
6457, 6471, 6472, 6478, 6479, 6504, 6507, 6509, 6532, 6536, 6538, 6552, 6558, 6559, 6621, 6624, 6662, 6667, 6669, 6670, 6785, 6787, 6788, 6789, 
6842, 6843, 6845, 6962, 6964, 6967, 6969, 7061, 7067, 7124, 7141, 7143, 7192, 7225, 7258, 7262, 7291, 7294, 7296, 7321, 7371, 7373, 7380, 7466, 
7490, 7504, 7511, 7516, 7517, 7518, 7519, 7524, 7571, 7572, 7573, 7578, 7579, 7674, 7680, 7681, 7781, 7784, 7803, 7805, 7807, 7809, 7815, 7915, 
7916, 7919, 7920, 7992, 7993, 7997, 7999, 8000, 8057, 8073, 8074, 8077, 8078, 8192, 8194, 8195, 8196, 8197, 8200, 8201, 8208, 8209, 8284, 8285, 
8289, 8301, 8302, 8304, 8306, 8308, 8350, 8402, 8404, 8405, 8410, 8433, 8472, 8474, 8485, 8553, 8604, 8610, 8621, 8625, 8626, 8627, 8628, 8630, 
8845, 8850, 8891, 8902, 8906, 8907, 8908, 8982, 8985, 8988, 9033, 9038, 9051, 9052, 9053, 9055, 9056, 9133, 9134, 9137, 9271, 9279, 9285, 9289, 
9297, 9299, 9306, 9315, 9332, 9333, 9336, 9338, 9346, 9421, 9422, 9424, 9428, 9429, 9430, 9464, 9466, 9546, 9547, 9556, 9557, 9558, 9559, 9687, 
9690, 9762, 9789, 9853, 9858, 9859, 9860, 9882, 9883, 9886, 9887, 9888, 9890, 10001, 10002, 10008, 10057, 10061, 10063, 10066, 10110, 10188, 10221, 
10321, 10323, 10325, 10422, 10424, 10425, 10426, 10427, 10431, 10433, 10434, 10436, 10438, 10452, 10643, 10647, 10686, 10712, 10714, 10715, 10718, 
10719, 10750, 10794, 10796, 10797, 10799, 10800, 10901, 10902, 10909, 10910, 10931, 10934, 10936, 10940, 11056, 11057, 11060, 11094, 11097, 11098, 
11100, 11267, 11363, 11364, 11365, 11366, 11367, 11370, 11392, 11394, 11398, 11399, 11403, 11404, 11405, 11408, 11413, 11457, 11491, 11494, 11497, 
11498, 11500, 11566, 11569, 11571, 11572, 11573, 11575, 11578, 11655, 11723, 11726, 11728, 11732, 11907, 11909, 11964, 11966, 11967, 11969, 11970, 
11973, 11976, 11977, 11980, 11986, 11988, 12050, 12085, 12111, 12112, 12113, 12114, 12115, 12116, 12117, 12118, 12119, 12120, 12143, 12144, 12145, 
12148, 12149, 12182, 12184, 12188, 12189, 12190, 12253, 12254, 12256, 12257, 12280, 12293, 12297, 12299, 12300, 12323, 12325, 12326, 12327, 12328, 
12330, 12341, 12344, 12345, 12347, 12348, 12351, 12353, 12358, 12359, 12385, 12387, 12388, 12452, 12453, 12456, 12457, 12458, 12466, 12468, 12469, 
12643, 12646, 12743, 12750, 12761, 12762, 12763, 12768, 12810, 12872, 12874, 12877, 12878, 12892, 12896, 12899, 12983, 13001, 13003, 13005, 13074, 
13080, 13161, 13162, 13181, 13182, 13188, 13189, 13190, 13266, 13268, 13291, 13294, 13295, 13296, 13300, 13392, 13397, 13398, 13399, 13400, 13401, 
13405, 13407, 13408, 13412, 13417, 13420, 13474, 13478, 13486, 13601, 13605, 13606, 13608, 13632, 13639, 13652, 13654, 13655, 13656, 13657, 13658, 
13659, 13660, 13711, 13712, 13713, 13715, 13716, 13717, 13718, 13720, 13751, 13752, 13753, 13754, 13755, 13758, 13780, 13802, 13803, 13809, 13832, 
13833, 13834, 13835, 13921, 13922, 13923, 13924, 13925, 13926, 13987, 13992, 13995, 13996, 13997, 13999, 14000, 14008, 14021, 14024, 14025, 14028, 
14029, 14183, 14189, 14251, 14258, 14281, 14283, 14285, 14286, 14287, 14333, 14335, 14340, 14372, 14373, 14375, 14379, 14380, 14382, 14385, 14401, 
14405, 14407, 14463, 14467, 14469, 14548, 14565, 14568, 14570, 14761, 14765, 14766, 14787, 14803, 14807, 14821, 14822, 14824, 15002, 15004, 15005, 
15007, 15008, 15036, 15072, 15077, 15078, 15110, 15117, 15120, 15131, 15134, 15137, 15138, 15139, 15196, 15198, 15203, 15207, 15222, 15231, 15235, 
15239, 15264, 15268, 15270, 15317, 15393, 15400, 15471, 15472, 15474, 15475, 15478, 15479, 15480, 15501, 15505, 15508, 15509, 15544, 15574, 15578, 
15579, 15621, 15625, 15626, 15627, 15630, 15747, 15748, 15750, 15796, 15852, 15853, 15857, 15859, 15860, 15943, 15946, 16015, 16017, 16022, 16024, 
16025, 16045, 16082, 16084, 16088, 16133, 16151, 16154, 16156, 16160, 16193, 16196, 16197, 16243, 16294, 16361, 16362, 16363, 16364, 16365, 16366, 
16367, 16368, 16369, 16370, 16438, 16439, 16567, 16601, 16602, 16603, 16604, 16605, 16606, 16607, 16608, 16609, 16610, 16615, 16616, 16618, 16651, 
16652, 16657, 16658, 16709, 16779, 16782, 16783, 16790, 16791, 16796, 16800, 16861, 16863, 16864, 16866, 16869, 16870, 16884, 16931, 16934, 16936, 
16944, 17004, 17041, 17048, 17221, 17224, 17225, 17226, 17227, 17228, 17229, 17230, 17330, 17349, 17411, 17536, 17540, 17547, 17590, 17611, 17613, 
17618, 17619, 17621, 17624, 17629, 17644, 17645, 17647, 17648, 17649, 17650, 17711, 17712, 17713, 17718, 17719, 17720, 17781, 17806, 17807, 17808, 
17809, 17830, 17841, 17871, 17875, 17891, 17892, 17893, 17900, 17901, 17902, 17907, 17909, 17944, 17991, 17998, 18037, 18123, 18124, 18125, 18126, 
18128, 18129, 18130, 18149, 18154, 18158, 18161, 18163, 18201, 18204, 18205, 18209, 18215, 18217, 18219, 18220, 18285, 18335, 18336, 18432, 18435, 
18440, 18446, 18472, 18473, 18475, 18476, 18477, 18481, 18484, 18485, 18487, 18489, 18490, 18553, 18555, 18557, 18560, 18621, 18622, 18623, 18624, 
18625, 18628, 18641, 18646, 18647, 18661, 18667, 18670, 18685, 18694, 18696, 18697, 18699, 18732, 18737, 18738, 18743, 18750, 18841, 18862, 18865, 
19052, 19059, 19091, 19093, 19094, 19098, 19100, 19141, 19142, 19143, 19148, 19150, 19189, 19234, 19236, 19243, 19291, 19292, 19293, 19295, 19300, 
19371, 19373, 19377, 19380, 19543, 19549, 19573, 19576, 19577, 19591, 19592, 19593, 19594, 19595, 19596, 19597, 19598, 19599, 19600, 19623, 19624, 
19625, 19626, 19628, 19630, 19771, 19774, 19776, 19777, 19783, 19785, 19786, 19787, 19793, 19800, 19852, 19859, 19901, 19904, 19907, 19909, 19941, 
19943, 19944, 19945, 19947, 19948, 19949, 19950, 19977, 19980, 19992, 19993, 19994, 19996, 19997, 19998, 20066, 20113, 20202, 20203, 20204, 20210, 
20329, 20332, 20338, 20504, 20509, 20510, 20532, 20535, 20538, 20539, 20544, 20547, 20552, 20555, 20556, 20560, 20680, 20690, 20782, 20784, 20787, 
20815, 20816, 20818, 20819, 20855, 20865, 20866, 20870, 20952, 20953, 20954, 20959, 21075, 21076, 21077, 21079, 21092, 21095, 21099, 21153, 21155, 
21159, 21292, 21293, 21295, 21296, 21298, 21302, 21342, 21343, 21345, 21346, 21550, 21601, 21602, 21603, 21604, 21605, 21606, 21607, 21608, 21609, 
21610, 21613, 21614, 21615, 21665, 21668, 21670, 21731, 21732, 21734, 21735, 21738, 21740, 21852, 21994, 21995, 21997, 22022, 22024, 22027, 22129, 
22130, 22440, 22500, 22593, 22597, 22664, 22666, 22668, 22670, 22691, 22694, 22695, 22732, 22736, 22740, 22773, 22919, 23110, 23157, 23186, 23247, 
23324, 23326, 23327, 23328, 23341, 23342, 23344, 23350, 23353, 23355, 23356, 23357, 23359, 23360, 23404, 23405, 23408, 23454, 23455, 23456, 23460, 
23475, 23492, 23494, 23495, 23577, 23578, 23580, 23608, 23652, 23702, 23704, 23723, 23724, 23727, 23730, 23791, 23796, 23813, 23815, 23816, 23847, 
23894, 23898, 23900, 23921, 23923, 23928, 23929, 23930, 23953, 23998, 24004, 24072, 24096, 24208, 24226, 24234, 24239, 24262, 24331, 24333, 24334, 
24338, 24340, 24724, 24725, 24727, 24741, 24842, 24861, 24974, 24982, 24986, 24987, 24988, 25026, 25099, 25111, 25114, 25157, 25160, 25207, 25210, 
25213, 25217, 25424, 25442, 25443, 25446, 25448, 25449, 25450, 25505, 25542, 25545, 25547, 25548, 25636, 25637, 25775, 25831, 25832, 25833, 25835, 
25836, 25838, 25901, 25922, 25925, 25936, 25981, 25982, 25984, 25985, 25986, 25987, 25988, 25990, 26012, 26013, 26016, 26018, 26019, 26108, 26111, 
26112, 26116, 26117, 26120, 26124, 26196, 26224, 26283, 26289, 26351, 26354, 26358, 26401, 26402, 26404, 26410, 26434, 26435, 26436, 26437, 26642, 
26643, 26646, 26647, 26649, 26781, 26881, 26943, 26947, 26950, 26968, 26969, 26984, 26986, 26987, 26989, 26990, 26992, 26999, 27173, 27292, 27294, 
27297, 27300, 27315, 27321, 27332, 27471, 27474, 27479, 27482, 27484, 27485, 27489, 27504, 27505, 27508, 27541, 27542, 27543, 27545, 27548, 27549, 
27551, 27554, 27560, 27571, 27574, 27576, 27578, 27585, 27741, 27760, 27880, 27901, 27902, 27904, 27909, 27910, 27974, 27996, 27999, 28000, 28052, 
28054, 28056, 28131, 28137, 28144, 28145, 28147, 28148, 28185, 28186, 28192, 28193, 28194, 28197, 28207, 28228, 28317, 28377, 28399, 28591, 28595, 
28596, 28598, 28599, 28650, 28751, 28752, 28753, 28754, 28760, 28827, 28855, 28867, 28891, 28893, 28895, 28900, 28957, 28961, 28962, 28965, 28966, 
28967, 28970, 28979, 28986, 29013, 29016, 29017, 29018, 29020, 29061, 29063, 29064, 29065, 29067, 29068, 29070, 29100, 29102, 29106, 29108, 29142, 
29143, 29144, 29147, 29159, 29294, 29296, 29297, 29299, 29374, 29378, 29391, 29395, 29398, 29441, 29443, 29446, 29447, 29448, 29480, 29543, 29544, 
29545, 29546, 29547, 29548, 29549, 29602, 29603, 29607, 29686, 29732, 29843, 29844, 29846, 29849, 29906, 29907, 29974, 29976, 29977, 29978, 29979, 
30001, 30002, 30003, 30008, 30009, 30010, 30014, 30015, 30016, 30018, 30019, 30071, 30072, 30073, 30076, 30077, 30078, 30171, 30173, 30175, 30181, 
30182, 30183, 30184, 30186, 30188, 30201, 30239, 30324, 30328, 30371, 30376, 30414, 30441, 30443, 30444, 30446, 30447, 30449, 30450, 30573, 30591, 
30593, 30599, 30600, 30647, 30653, 30660, 30717, 30733, 30734, 30743, 30744, 30803, 30805, 30807, 30851, 30858, 30881, 30882, 30886, 30889, 30890, 
30942, 30961, 30963, 30968, 31025, 31029, 31131, 31132, 31133, 31134, 31136, 31137, 31138, 31139, 31140, 31248, 31271, 31277, 31278, 31292, 31294, 
31298, 31299, 31300, 31302, 31309, 31326, 31328, 31334, 31383, 31385, 31420, 31513, 31514, 31515, 31518, 31519, 31520, 31602, 31651, 31657, 31659, 
31660, 31682, 31685, 31687, 31689, 31864, 31873, 31875, 31877, 31878, 31880, 31882, 31951, 31956, 31957, 31958, 31960, 32015, 32063, 32064, 32066, 
32069, 32088, 32089, 32090, 32181, 32211, 32213, 32215, 32216, 32217, 32218, 32251, 32252, 32256, 32258, 32259, 32260, 32271, 32272, 32273, 32274, 
32276, 32278, 32279, 32284, 32286, 32334, 32340, 32393, 32434, 32441, 32442, 32443, 32448, 32473, 32474, 32476, 32480, 32523, 32524, 32526, 32530, 
32537, 32540, 32674, 32676, 32680, 32696, 32698, 32700, 32701, 32707, 32708, 32710, 32748, 32908, 32909, 32961, 32965, 32966, 32969, 32970, 33041, 
33043, 33047, 33048, 33235, 33237, 33241, 33261, 33268, 33269, 33283, 33326, 33409, 33472, 33475, 33478, 33637, 33652, 33655, 33656, 33659, 33781, 
33782, 33784, 33785, 33801, 33802, 33804, 33805, 33806, 33809, 33877, 33879, 33931, 33932, 33933, 33936, 33938, 33940, 33976, 34063, 34115, 34287, 
34311, 34317, 34320, 34353, 34355, 34357, 34359, 34360, 34566, 34613, 34615, 34620, 34722, 34729, 34751, 34753, 34755, 34756, 34760, 34762, 34833, 
34835, 34841, 34848, 34903, 34916, 34984, 34985, 35133, 35165, 35169, 35181, 35182, 35183, 35184, 35185, 35186, 35187, 35188, 35189, 35190, 35202, 
35206, 35207, 35268, 35336, 35342, 35376, 35379, 35382, 35383, 35384, 35387, 35389, 35390, 35405, 35516, 35534, 35538, 35540, 35662, 35663, 35666, 
35668, 35669, 35670, 35691, 35693, 35694, 35711, 35712, 35714, 35717, 35718, 35721, 35724, 35726, 35727, 35728, 35730, 35798, 35806, 35833, 35904, 
35938, 35969, 36025, 36029, 36143, 36163, 36169, 36175, 36204, 36229, 36465, 36480, 36491, 36498, 36532, 36571, 36572, 36574, 36575, 36578, 36579, 
36582, 36583, 36584, 36585, 36586, 36587, 36622, 36623, 36624, 36625, 36627, 36629, 36643, 36647, 36648, 36649, 36650, 36681, 36775, 36779, 36844, 
36845, 36848, 36849, 36858, 36861, 36864, 36868, 36870, 36875, 36925, 37001, 37004, 37007, 37008, 37010, 37142, 37148, 37154, 37158, 37159, 37161, 
37163, 37169, 37352, 37354, 37360, 37375, 37377, 37378, 37381, 37426, 37427, 37441, 37442, 37445, 37447, 37449, 37450, 37461, 37463, 37466, 37523, 
37528, 37543, 37546, 37548, 37549, 37631, 37632, 37633, 37636, 37637, 37638, 37639, 37640, 37649, 37691, 37692, 37699, 37701, 37702, 37703, 37704, 
37705, 37710, 37721, 37728, 37730, 37742, 37746, 37748, 37827, 37873, 37952, 37953, 37954, 37955, 37956, 37957, 37958, 37960, 37979, 37991, 38000, 
38032, 38033, 38035, 38036, 38037, 38038, 38085, 38086, 38087, 38088, 38121, 38122, 38123, 38125, 38129, 38130, 38141, 38146, 38147, 38148, 38149, 
38281, 38282, 38284, 38286, 38287, 38301, 38303, 38305, 38307, 38309, 38372, 38373, 38375, 38377, 38379, 38380, 38394, 38395, 38396, 38431, 38433, 
38435, 38439, 38476, 38544, 38561, 38562, 38563, 38566, 38570, 38575, 38578, 38579, 38580, 38586, 38587, 38589, 38590, 38681, 38683, 38684, 38688, 
38689, 38690, 38713, 38811, 38812, 38813, 38818, 38819, 38971, 38977, 38978, 38980, 39001, 39003, 39004, 39010, 39093, 39094, 39161, 39162, 39163, 
39165, 39166, 39167, 39168, 39169, 39170, 39244, 39262, 39266, 39283, 39285, 39301, 39307, 39308, 39310, 39315, 39317, 39319, 39364, 39414, 39419, 
39441, 39444, 39449, 39450, 39474, 39481, 39482, 39483, 39484, 39485, 39486, 39488, 39490, 39502, 39510, 39532, 39533, 39537, 39538, 39541, 39545, 
39546, 39549, 39550, 39567, 39569, 39570, 39658, 39659, 39776, 39857, 39904, 39973, 39975, 39978, 40048, 40053, 40081, 40082, 40084, 40085, 40088, 
40120, 40140, 40156, 40159, 40197, 40204, 40223, 40226, 40237, 40273, 40275, 40276, 40278, 40291, 40363, 40364, 40367, 40370, 40378, 40379, 40380, 
40542, 40545, 40546, 40548, 40549, 40586, 40590, 40661, 40662, 40663, 40666, 40667, 40670, 40695, 40699, 40755, 40756, 40757, 40759, 40794, 40795, 
40798, 40799, 40804, 40821, 40826, 40841, 40846, 40847, 40894, 40897, 40898, 40908, 40911, 40914, 40961, 40962, 40964, 40967, 40968, 40977, 41016, 
41128, 41188, 41238, 41324, 41327, 41328, 41342, 41343, 41347, 41348, 41353, 41448, 41474, 41476, 41477, 41479, 41512, 41516, 41517, 41518, 41520, 
41521, 41523, 41524, 41525, 41527, 41530, 41624, 41626, 41627, 41628, 41629, 41634, 41635, 41639, 41726, 41751, 41758, 41765, 41766, 41767, 41768, 
41770, 41834, 41835, 41837, 41981, 41982, 41983, 41985, 42013, 42017, 42052, 42084, 42143, 42144, 42145, 42146, 42147, 42148, 42149, 42172, 42180, 
42251, 42257, 42260, 42292, 42295, 42326, 42329, 42371, 42372, 42373, 42374, 42375, 42376, 42377, 42378, 42379, 42380, 42402, 42404, 42406, 42409, 
42410, 42471, 42475, 42477, 42478, 42479, 42480, 42484, 42485, 42486, 42487, 42488, 42489, 42490, 42492, 42494, 42495, 42496, 42498, 42499, 42556, 
42558, 42627, 42628, 42629, 42653, 42654, 42657, 42807, 42934, 42991, 42994, 42996, 42997, 43016, 43019, 43025, 43032, 43037, 43038, 43040, 43042, 
43043, 43044, 43047, 43048, 43050, 43145, 43150, 43323, 43337, 43361, 43363, 43367, 43370, 43487, 43495, 43496, 43498, 43499, 43500, 43523, 43524, 
43525, 43526, 43527, 43530, 43599, 43616, 43618, 43619, 43621, 43625, 43629, 43774, 43775, 43777, 43781, 43846, 43922, 43924, 43927, 43928, 43930, 
43962, 44011, 44014, 44017, 44020, 44077, 44083, 44090, 44092, 44161, 44162, 44164, 44165, 44167, 44169, 44382, 44388, 44389, 44516, 44603, 44605, 
44606, 44607, 44612, 44614, 44616, 44631, 44647, 44651, 44657, 44772, 44781, 44788, 44790, 44802, 44806, 44807, 44808, 44967, 45052, 45053, 45096, 
45115, 45140, 45151, 45358, 45359, 45394, 45397, 45471, 45472, 45474, 45475, 45477, 45478, 45479, 45480, 45501, 45502, 45503, 45504, 45505, 45506, 
45507, 45508, 45509, 45510, 45582, 45583, 45584, 45645, 45649, 45661, 45662, 45668, 45669, 45716, 45805, 45808, 45809, 45911, 45913, 45914, 45915, 
45917, 45918, 45920, 45926, 45967, 45993, 46044, 46045, 46046, 46050, 46061, 46062, 46069, 46070, 46101, 46102, 46103, 46104, 46108, 46109, 46125, 
46126, 46129, 46130, 46275, 46278, 46280, 46301, 46305, 46306, 46307, 46388, 46419, 46559, 46601, 46604, 46605, 46607, 46610, 46675, 46676, 46677, 
46695, 46698, 46700, 46701, 46702, 46703, 46704, 46707, 46709, 46710, 46821, 46823, 46824, 46826, 46830, 46833, 46867, 46868, 46991, 46994, 46995, 
46996, 46997, 46998, 46999, 47000, 47009, 47184, 47192, 47194, 47274, 47275, 47279, 47282, 47297, 47444, 47490, 47492, 47493, 47502, 47508, 47509, 
47510, 47588, 47589, 47681, 47682, 47687, 47688, 47690, 47691, 47692, 47694, 47695, 47697, 47698, 47699, 47700, 47753, 47754, 47756, 47757, 47758, 
47802, 47803, 47806, 47820, 47910, 47925, 47927, 47929, 48011, 48012, 48013, 48018, 48042, 48048, 48049, 48064, 48066, 48068, 48094, 48095, 48096, 
48097, 48158, 48171, 48172, 48173, 48176, 48177, 48178, 48182, 48253, 48255, 48256, 48258, 48347, 48348, 48349, 48372, 48374, 48375, 48376, 48401, 
48403, 48405, 48409, 48410, 48454, 48459, 48532, 48540, 48557, 48692, 48707, 48708, 48760, 48772, 48864, 48866, 48867, 48870, 48871, 48875, 48876, 
48877, 48878, 48921, 48923, 48926, 48927, 48928, 48932, 48939, 48940, 49131, 49132, 49134, 49135, 49140, 49214, 49220, 49221, 49222, 49227, 49230, 
49350, 49501, 49502, 49503, 49504, 49507, 49525, 49526, 49527, 49530, 49555, 49556, 49557, 49612, 49762, 49771, 49778, 49782, 49785, 49796, 49797, 
49798, 49799, 49853, 49860, 50134, 50139, 50176, 50201, 50206, 50209, 50263, 50286, 50289, 50301, 50303, 50304, 50307};

void fill_hot_dead_map_eta_bin_centers(int runnumber, float minus_z, float plus_z) {

	// hot dead maps 
    vector<int> *emcal_hot_dead_ieta = nullptr;
    vector<int> *emcal_hot_dead_iphi = nullptr;
    vector<int> *ihcal_hot_dead_ieta = nullptr;
    vector<int> *ihcal_hot_dead_iphi = nullptr;
    vector<int> *ohcal_hot_dead_ieta = nullptr;
    vector<int> *ohcal_hot_dead_iphi = nullptr;
    vector<float> *ihcal_eta_bins = nullptr;
    vector<float> *ohcal_eta_bins = nullptr;

	//TFile *hotdeadfile = new TFile(TString::Format("run%d_hotdeadmap_z_%d_%d_p011.root", runnumber, int(floor(minus_z)), int(floor(plus_z))), "READ"); // edited for new data
	TFile *hotdeadfile = new TFile(TString::Format("/sphenix/user/egm2153/calib_study/detdeta/analysis/runs23727_23746/run%d_hotdeadmap_z_%d_%d_new_status.root", runnumber, int(floor(minus_z)), int(floor(plus_z))), "READ");
	TTree *hotdeadtree = dynamic_cast<TTree*>(hotdeadfile->Get("T"));
	hotdeadtree->SetBranchAddress("emcal_hot_dead_ieta", &emcal_hot_dead_ieta);
    hotdeadtree->SetBranchAddress("emcal_hot_dead_iphi", &emcal_hot_dead_iphi);
    hotdeadtree->SetBranchAddress("ihcal_hot_dead_ieta", &ihcal_hot_dead_ieta);
    hotdeadtree->SetBranchAddress("ihcal_hot_dead_iphi", &ihcal_hot_dead_iphi);
    hotdeadtree->SetBranchAddress("ohcal_hot_dead_ieta", &ohcal_hot_dead_ieta);
    hotdeadtree->SetBranchAddress("ohcal_hot_dead_iphi", &ohcal_hot_dead_iphi);

    // edited 3.25.24 to use another run in run group 
    hotdeadtree->GetEntry(0);

	//std::cout << "emcal hot dead map input size " << emcal_hot_dead_ieta->size() << std::endl;
	for (int i = 0; i < emcal_hot_dead_ieta->size(); i++) {
		tuple<int, int> new_hot_tower = make_tuple((*emcal_hot_dead_ieta)[i], (*emcal_hot_dead_iphi)[i]); 
        emcal_hot_dead_map.insert(new_hot_tower);
	}
	//std::cout << "emcal hot dead map set size " << emcal_hot_dead_map.size() << std::endl;
	
	//std::cout << "ihcal hot dead map input size " << ihcal_hot_dead_ieta->size() << std::endl;
	for (int i = 0; i < ihcal_hot_dead_ieta->size(); i++) {
		tuple<int, int> new_hot_tower = make_tuple((*ihcal_hot_dead_ieta)[i], (*ihcal_hot_dead_iphi)[i]); 
        ihcal_hot_dead_map.insert(new_hot_tower);
	}
	//std::cout << "ihcal hot dead map set size " << ihcal_hot_dead_map.size() << std::endl;

	//std::cout << "ohcal hot dead map input size " << ohcal_hot_dead_ieta->size() << std::endl;
	for (int i = 0; i < ohcal_hot_dead_ieta->size(); i++) {
		tuple<int, int> new_hot_tower = make_tuple((*ohcal_hot_dead_ieta)[i], (*ohcal_hot_dead_iphi)[i]); 
        ohcal_hot_dead_map.insert(new_hot_tower);
	}
	//std::cout << "ohcal hot dead map set size " << ohcal_hot_dead_map.size() << std::endl;
	hotdeadfile->Close();

	//TFile *etabinfile = new TFile(TString::Format("run23727_hotdeadmap_z_%d_%d_p011.root", int(floor(minus_z)), int(floor(plus_z))), "READ"); // edited for new data
	TFile *etabinfile = new TFile(TString::Format("/sphenix/user/egm2153/calib_study/detdeta/analysis/runs23727_23746/run23727_hotdeadmap_z_%d_%d_new_status.root", int(floor(minus_z)), int(floor(plus_z))), "READ"); 
	TTree *etabintree = dynamic_cast<TTree*>(etabinfile->Get("T"));
	
	etabintree->SetBranchAddress("ihcal_eta_bin_centers", &ihcal_eta_bins);
    etabintree->SetBranchAddress("ohcal_eta_bin_centers", &ohcal_eta_bins);
	etabintree->GetEntry(0);
	for (int i = 0; i < ihcal_eta_bins->size(); i++) {
		ihcal_eta_bin_centers[i] = (*ihcal_eta_bins)[i];
		emcal_eta_bin_centers[i] = (*ihcal_eta_bins)[i];
		calo_eta_bin_centers[i] = (*ihcal_eta_bins)[i];
	}
	for (int i = 0; i < ohcal_eta_bins->size(); i++) {
		ohcal_eta_bin_centers[i] = (*ohcal_eta_bins)[i];
	}

	etabinfile->Close();

}

float mc_centrality(int centbin) {
	int centrality_index = centbin/5;
	return centrality_bin[19-centrality_index];
}

void fill_zvertex_centrality(int dataormc, int runnumber, const char* generator) {
	
	TFile *zvertexfile;
	float vz_reweight[200]; 
	float mc_cent[20];
	float data_cent[20];
	if (dataormc) {
		zvertexfile = new TFile(TString::Format("/sphenix/user/egm2153/calib_study/detdeta/analysis/runs23727_23746/dETdeta_vertex_reweight_run%d_reweight_hijing_new_status.root", runnumber), "READ");
		TTree *vertextree = dynamic_cast<TTree*>(zvertexfile->Get("T"));
		vertextree->SetBranchAddress("vertex_reweight", vz_reweight);
		vertextree->SetBranchAddress("mc_centrality", mc_cent);
		vertextree->GetEntry(0);
		for (int i = 0; i < 200; i++) {
			vertex_reweight[i] = vz_reweight[i];
		}
		centrality_bin.assign(mc_cent, mc_cent+20);
	} else { // edited 3.25.24 from epos to reweight_hijing
		zvertexfile = new TFile(TString::Format("/sphenix/user/egm2153/calib_study/detdeta/analysis/runs23727_23746/dETdeta_vertex_reweight_run%d_reweight_hijing_new_status.root",runnumber), "READ");
		TTree *vertextree = dynamic_cast<TTree*>(zvertexfile->Get("T"));
		vertextree->SetBranchAddress("data_centrality", data_cent);
		vertextree->GetEntry(0);
		centrality_bin.assign(data_cent, data_cent+20);
	}
	zvertexfile->Close();

}

void fill_zs_cross_calib() {

	float emcal_zs_calib[96][256];
	float ihcal_zs_calib[24][64];
	float ohcal_zs_calib[24][64];
	TFile *zscalibfile = new TFile("/sphenix/user/egm2153/calib_study/zs_testing/zs_testing_run23727_p008_new.root", "READ");
	TTree *zscalibtree = dynamic_cast<TTree*>(zscalibfile->Get("zs_calib_tree"));
	zscalibtree->SetBranchAddress("emcal_zs_calib", emcal_zs_calib);
	zscalibtree->SetBranchAddress("ihcal_zs_calib", ihcal_zs_calib);
	zscalibtree->SetBranchAddress("ohcal_zs_calib", ohcal_zs_calib);
	zscalibtree->GetEntry(0);
	for (int i = 0; i < 96; i++) {
		for (int j = 0; j < 256; j++) {
			em_zs_calib[i][j] = emcal_zs_calib[i][j];
			if (i < 24 && j < 64) {
				ih_zs_calib[i][j] = ihcal_zs_calib[i][j];
				oh_zs_calib[i][j] = ohcal_zs_calib[i][j];
			}
		}
	}
	zscalibfile->Close();
}

void dETdeta_analysis(int runnumber = 23727, const char* generator = "", float minus_z = -2, float plus_z = 2, int dataormc = 0, int reweighting = 1, int central = 0, int min_cent = 0, int max_cent = 5, int zs = 0, int zs_value = 10, int time = 0, const char* opt_tag = "") {

	string filename = "dETdeta_analysis_";
	filename += to_string(runnumber);
	string opttag = opt_tag;
	if (strcmp(opt_tag,"")) { filename += "_" + opttag; }
	string zstag; 
	if (zs == 0) { zstag = "nozs"; }
	if (zs == 1) { zstag = "zs_" + to_string(zs_value) + "ADC"; }
	if (zs == 2) { zstag = "zs_ecut" + to_string(zs_value) + "ADC"; }
	string timetag = (time?"emcal_timecut_":"");
  	string dattag = (dataormc?"mc":"data");
  	string weighttag = (reweighting?"reweight":"noweight");
  	string centtag = (central?to_string(min_cent)+"-"+to_string(max_cent):"0-90"); 
  	string gentag = generator;
  	if (!strcmp(generator,"")) filename += "_" + zstag + "_" + timetag + dattag + "_" + weighttag + "_" + centtag + ".root";
  	else filename += "_" + zstag + "_" + timetag + dattag + "_" + weighttag + "_" + centtag + "_" + gentag + ".root";

	fill_hot_dead_map_eta_bin_centers(runnumber, minus_z, plus_z);
	fill_zvertex_centrality(dataormc, runnumber, generator);
	fill_zs_cross_calib();
	assert(centrality_bin.size() == 20);

	std::cout << "centrality_bins" << std::endl;
	for (int i = 0; i < 20; i++) {
		std::cout << centrality_bin[i] << " ";
	}
	std::cout << std::endl;

	TFile *out = new TFile(filename.c_str(),"RECREATE");
	TH1F* h_vz = new TH1F("h_vz","",400, -100, 100);
	TH1F* h_vz_reweight = new TH1F("h_vz_reweight","",400, -100, 100);
	TH1F* h_mbd = new TH1F("h_mbd","",250000,0,250000);

	TH1F* h_event_truth_energy = new TH1F("h_event_truth_energy","",10000,0,10000);
	TH1F* hetdeta = new TH1F("hetdeta","",240,-12,12);
	TH1F* hetdeta_zoom = new TH1F("hetdeta_zoom","",220,-1.1,1.1);
	
	TH2F* h_2D_ihcal_calib = new TH2F("h_2D_ihcal_calib","",24,0.,24.,64,0.,64.);
  	TH2F* h_2D_ohcal_calib = new TH2F("h_2D_ohcal_calib","",24,0.,24.,64,0.,64.);
  	TH2F* h_2D_emcal_calib = new TH2F("h_2D_emcal_calib","",96,0.,96.,256,0.,256.);

  	TH2F* h_2D_ihcal_calibT = new TH2F("h_2D_ihcal_calibT","",24,0.,24.,64,0.,64.);
  	TH2F* h_2D_ohcal_calibT = new TH2F("h_2D_ohcal_calibT","",24,0.,24.,64,0.,64.);
  	TH2F* h_2D_emcal_calibT = new TH2F("h_2D_emcal_calibT","",96,0.,96.,256,0.,256.);
	
	TH1F* h_event_energy = new TH1F("h_event_energy","", 6000,-1000,5000);
	TH1F* h_event_hcal_energy = new TH1F("h_event_hcal_energy","", 6000,-1000,5000);
	TH1F* h_event_emcal_energy = new TH1F("h_event_emcal_energy","", 6000,-1000,5000);
	TH1F* h_event_ihcal_energy = new TH1F("h_event_ihcal_energy","", 2000,-1000,1000);
	TH1F* h_event_ohcal_energy = new TH1F("h_event_ohcal_energy","", 2000,-1000,1000);

	TH1F* h_emcal = new TH1F("h_emcal","",1000,0,10);
	TH1F* h_ihcal = new TH1F("h_ihcal","",1000,0,10);
	TH1F* h_ohcal = new TH1F("h_ohcal","",1000,0,10);

	TH2F* h_em_zero_zscross = new TH2F("h_em_zero_zscross","",96,0,96,256,0,256);
	TH2F* h_ih_zero_zscross = new TH2F("h_ih_zero_zscross","",24,0,24,64,0,64);
	TH2F* h_oh_zero_zscross = new TH2F("h_oh_zero_zscross","",24,0,24,64,0,64);


	
	int emcal_num_bins = 24;
	double emcal_bin_edges[25];
    if (!dataormc || reweighting) {
    	emcal_bin_edges[0] = emcal_eta_bin_centers[0] - 0.5 * (emcal_eta_bin_centers[1] - emcal_eta_bin_centers[0]);
	    for (int i = 1; i < emcal_num_bins; ++i) { emcal_bin_edges[i] = (emcal_eta_bin_centers[i] + emcal_eta_bin_centers[i - 1]) / 2.0; }
	    emcal_bin_edges[emcal_num_bins] = emcal_eta_bin_centers[emcal_num_bins - 1] + 0.5 * (emcal_eta_bin_centers[emcal_num_bins - 1] - emcal_eta_bin_centers[emcal_num_bins - 2]);
    } else {
    	emcal_bin_edges[0] = eta_bin_centers[0] - 0.5 * (eta_bin_centers[1] - eta_bin_centers[0]);
    	for (int i = 1; i < emcal_num_bins; ++i) { emcal_bin_edges[i] = (eta_bin_centers[i] + eta_bin_centers[i - 1]) / 2.0; }
    	emcal_bin_edges[emcal_num_bins] = eta_bin_centers[emcal_num_bins - 1] + 0.5 * (eta_bin_centers[emcal_num_bins - 1] - eta_bin_centers[emcal_num_bins - 2]);
    }
    TH1F* h_eT_emcal = new TH1F("h_eT_emcal","",emcal_num_bins,emcal_bin_edges);
	TH1F* hetdeta_emcalbin = new TH1F("hetdeta_emcalbin","",emcal_num_bins, emcal_bin_edges);
	
	int ihcal_num_bins = 24;
    double ihcal_bin_edges[25];
    if (!dataormc || reweighting) {
    	ihcal_bin_edges[0] = ihcal_eta_bin_centers[0] - 0.5 * (ihcal_eta_bin_centers[1] - ihcal_eta_bin_centers[0]);
	    for (int i = 1; i < ihcal_num_bins; ++i) { ihcal_bin_edges[i] = (ihcal_eta_bin_centers[i] + ihcal_eta_bin_centers[i - 1]) / 2.0; }
	    ihcal_bin_edges[ihcal_num_bins] = ihcal_eta_bin_centers[ihcal_num_bins - 1] + 0.5 * (ihcal_eta_bin_centers[ihcal_num_bins - 1] - ihcal_eta_bin_centers[ihcal_num_bins - 2]);
    } else {
    	ihcal_bin_edges[0] = eta_bin_centers[0] - 0.5 * (eta_bin_centers[1] - eta_bin_centers[0]);
    	for (int i = 1; i < ihcal_num_bins; ++i) { ihcal_bin_edges[i] = (eta_bin_centers[i] + eta_bin_centers[i - 1]) / 2.0; }
    	ihcal_bin_edges[ihcal_num_bins] = eta_bin_centers[ihcal_num_bins - 1] + 0.5 * (eta_bin_centers[ihcal_num_bins - 1] - eta_bin_centers[ihcal_num_bins - 2]);
    }
    TH1F* h_eT_ihcal = new TH1F("h_eT_ihcal","",ihcal_num_bins,ihcal_bin_edges);
	TH1F* hetdeta_ihcalbin = new TH1F("hetdeta_ihcalbin","",ihcal_num_bins, ihcal_bin_edges);
	
	int ohcal_num_bins = 24;
    double ohcal_bin_edges[25];
    if (!dataormc || reweighting) {
	    ohcal_bin_edges[0] = ohcal_eta_bin_centers[0] - 0.5 * (ohcal_eta_bin_centers[1] - ohcal_eta_bin_centers[0]);
	    for (int i = 1; i < ohcal_num_bins; ++i) { ohcal_bin_edges[i] = (ohcal_eta_bin_centers[i] + ohcal_eta_bin_centers[i - 1]) / 2.0; }
	    ohcal_bin_edges[ohcal_num_bins] = ohcal_eta_bin_centers[ohcal_num_bins - 1] + 0.5 * (ohcal_eta_bin_centers[ohcal_num_bins - 1] - ohcal_eta_bin_centers[ohcal_num_bins - 2]);
	} else {
		ohcal_bin_edges[0] = eta_bin_centers[0] - 0.5 * (eta_bin_centers[1] - eta_bin_centers[0]);
	    for (int i = 1; i < ohcal_num_bins; ++i) { ohcal_bin_edges[i] = (eta_bin_centers[i] + eta_bin_centers[i - 1]) / 2.0; }
	    ohcal_bin_edges[ohcal_num_bins] = eta_bin_centers[ohcal_num_bins - 1] + 0.5 * (eta_bin_centers[ohcal_num_bins - 1] - eta_bin_centers[ohcal_num_bins - 2]);
	}
	TH1F* h_eT_ohcal = new TH1F("h_eT_ohcal","",ohcal_num_bins,ohcal_bin_edges);
	TH1F* hetdeta_ohcalbin = new TH1F("hetdeta_ohcalbin","",ohcal_num_bins, ohcal_bin_edges);
	
	int calo_num_bins = 24;
    double calo_bin_edges[25];
    if (!dataormc || reweighting) {
    	calo_bin_edges[0] = calo_eta_bin_centers[0] - 0.5 * (calo_eta_bin_centers[1] - calo_eta_bin_centers[0]);
	    for (int i = 1; i < calo_num_bins; ++i) { calo_bin_edges[i] = (calo_eta_bin_centers[i] + calo_eta_bin_centers[i - 1]) / 2.0; }
	    calo_bin_edges[calo_num_bins] = calo_eta_bin_centers[calo_num_bins - 1] + 0.5 * (calo_eta_bin_centers[calo_num_bins - 1] - calo_eta_bin_centers[calo_num_bins - 2]);
    } else {
    	calo_bin_edges[0] = eta_bin_centers[0] - 0.5 * (eta_bin_centers[1] - eta_bin_centers[0]);
    	for (int i = 1; i < calo_num_bins; ++i) { calo_bin_edges[i] = (eta_bin_centers[i] + eta_bin_centers[i - 1]) / 2.0; }
    	calo_bin_edges[calo_num_bins] = eta_bin_centers[calo_num_bins - 1] + 0.5 * (eta_bin_centers[calo_num_bins - 1] - eta_bin_centers[calo_num_bins - 2]);
    }
    TH1F* h_eT_calo = new TH1F("h_eT_calo","",calo_num_bins,calo_bin_edges);
	TH1F* hetdeta_calobin = new TH1F("hetdeta_calobin","",calo_num_bins, calo_bin_edges);
	
    TChain chain("ttree");

    if (dataormc && !strcmp(generator, "ftfp_bert")) {
	    TString wildcardPath = "/sphenix/user/egm2153/calib_study/detdeta/eposrun/ftfp_bert_condor/epos_sim_output.root";
	    chain.Add(wildcardPath);
    } else if (dataormc && !strcmp(generator, "ftfp_bert_hp")) { 
	    TString wildcardPath = "/sphenix/user/egm2153/calib_study/detdeta/eposrun/ftfp_bert_hp_condor/epos_sim_output.root";
	    chain.Add(wildcardPath);
    } else if (dataormc && !strcmp(generator, "qgsp_bert_hp")) {
	    TString wildcardPath = "/sphenix/user/egm2153/calib_study/detdeta/eposrun/qgsp_bert_hp_condor/epos_sim_output.root";
	    chain.Add(wildcardPath);
    } else {
    	std::cout << "generator/data not found" << std::endl;
    	return;
    }

    int m_simtwrmult_cemc;
    float m_simtwr_cemc_e[cemcSize];
    float m_simtwr_cemc_zs_e[cemcSize];
    int m_simtwr_cemc_ieta[cemcSize];
    int m_simtwr_cemc_iphi[cemcSize];
    float m_simtwr_cemc_eta[cemcSize];
    int m_simtwr_cemc_adc[cemcSize];
    int m_simtwr_cemc_zs_adc[cemcSize];
    float m_simtwr_cemc_time[cemcSize];
    int m_simtwrmult_ihcal;
    float m_simtwr_ihcal_e[ihcalSize];
    float m_simtwr_ihcal_zs_e[ihcalSize];
    int m_simtwr_ihcal_ieta[ihcalSize];
    int m_simtwr_ihcal_iphi[ihcalSize];
    float m_simtwr_ihcal_eta[ihcalSize];
    int m_simtwr_ihcal_adc[ihcalSize];
    int m_simtwr_ihcal_zs_adc[ihcalSize];
    float m_simtwr_ihcal_time[ihcalSize];
    int m_simtwrmult_ohcal;
    float m_simtwr_ohcal_e[ohcalSize];
    float m_simtwr_ohcal_zs_e[ohcalSize];
    int m_simtwr_ohcal_ieta[ohcalSize];
    int m_simtwr_ohcal_iphi[ohcalSize];
    float m_simtwr_ohcal_eta[ohcalSize];
    int m_simtwr_ohcal_adc[ohcalSize];
    int m_simtwr_ohcal_zs_adc[ohcalSize];
    float m_simtwr_ohcal_time[ohcalSize];
    int m_sectormb;
    float m_mbenergy[mbdSize];
    int m_g4;
    float m_g4_e[g4Size];
    float m_g4_eta[g4Size];
    float m_g4_pt[g4Size];
    float m_g4_pz[g4Size];
    float m_vtx[vtxSize];
    bool m_isMinBias;
    int m_centbin;

     // Set branch addresses
    chain.SetBranchAddress("sectorem", &m_simtwrmult_cemc);
    chain.SetBranchAddress("emcalen", m_simtwr_cemc_e);
    chain.SetBranchAddress("emcaletabin", m_simtwr_cemc_ieta);
    chain.SetBranchAddress("emcalphibin", m_simtwr_cemc_iphi);
    chain.SetBranchAddress("emetacor", m_simtwr_cemc_eta);
    if (!dataormc) {
    	chain.SetBranchAddress("emcalzs", m_simtwr_cemc_zs_e);
	    chain.SetBranchAddress("emcaladc", m_simtwr_cemc_adc);
	    chain.SetBranchAddress("emcalzsadc", m_simtwr_cemc_zs_adc);
	    chain.SetBranchAddress("emcalt", m_simtwr_cemc_time);
	}
    chain.SetBranchAddress("sectorih", &m_simtwrmult_ihcal);
    chain.SetBranchAddress("ihcalen", m_simtwr_ihcal_e);
    chain.SetBranchAddress("ihcaletabin", m_simtwr_ihcal_ieta);
    chain.SetBranchAddress("ihcalphibin", m_simtwr_ihcal_iphi);
    chain.SetBranchAddress("ihetacor", m_simtwr_ihcal_eta);
    if (!dataormc) {
    	chain.SetBranchAddress("ihcalzs", m_simtwr_ihcal_zs_e);
	    chain.SetBranchAddress("ihcaladc", m_simtwr_ihcal_adc);
	    chain.SetBranchAddress("ihcalzsadc", m_simtwr_ihcal_zs_adc);
	    chain.SetBranchAddress("ihcalt", m_simtwr_ihcal_time);
	}
    chain.SetBranchAddress("sectoroh", &m_simtwrmult_ohcal);
    chain.SetBranchAddress("ohcalen", m_simtwr_ohcal_e);
    chain.SetBranchAddress("ohcaletabin", m_simtwr_ohcal_ieta);
    chain.SetBranchAddress("ohcalphibin", m_simtwr_ohcal_iphi);
    chain.SetBranchAddress("ohetacor", m_simtwr_ohcal_eta);
    if (!dataormc) {
    	chain.SetBranchAddress("ohcalzs", m_simtwr_ohcal_zs_e);
	    chain.SetBranchAddress("ohcaladc", m_simtwr_ohcal_adc);
	    chain.SetBranchAddress("ohcalzsadc", m_simtwr_ohcal_zs_adc);
	    chain.SetBranchAddress("ohcalt", m_simtwr_ohcal_time);
	}
    chain.SetBranchAddress("sectormb", &m_sectormb);
    chain.SetBranchAddress("mbenrgy", &m_mbenergy);
    if (!dataormc) {
    	chain.SetBranchAddress("isMinBias", &m_isMinBias);
	    chain.SetBranchAddress("centbin", &m_centbin);
    }
    if (dataormc) {
    	chain.SetBranchAddress("truthpar_n", &m_g4);
	    chain.SetBranchAddress("truthpar_e", m_g4_e);
	    chain.SetBranchAddress("truthpar_eta", m_g4_eta);
	    chain.SetBranchAddress("truthpar_pt", m_g4_pt);
	    chain.SetBranchAddress("truthpar_pz", m_g4_pz);
    }
    chain.SetBranchAddress("track_vtx", m_vtx);

	int eventnumber = 0;
	float totalweights = 0.0;

    Long64_t nEntries = chain.GetEntries();
    std::cout << nEntries << std::endl;
    for (Long64_t entry = 0; entry < nEntries; ++entry) {
    //for (Long64_t entry = 0; entry < 40000; ++entry) {
        chain.GetEntry(entry);
    	if (eventnumber % 1000 == 0) cout << "event " << eventnumber << endl;

		float emcale = 0;
		float ihcale = 0;
		float ohcale = 0;
		float totale = 0;
		float truthe = 0;

		float E_emcal[96] = {0};
		float E_ihcal[24] = {0};
  		float E_ohcal[24] = {0};

  		float vz_weight = 1.0;

  		eventnumber++;

  		if (!strcmp(generator, "ftfp_bert_hp") || !strcmp(generator, "qgsp_bert_hp")) {
  			auto sit = skip_event.find(eventnumber);
		    if (sit != skip_event.end()) { continue; }
  		}

  		// require that simulation could reconstruct a vertex for the event
  		//if(isnan(m_vtx[2])) { continue; }
  		//if (m_vtx[2] < minus_z || m_vtx[2] > plus_z) { continue; }
  		if (!dataormc && !m_isMinBias) { continue; }

  		float totalcharge = 0.0;
  		for (int i = 0; i < m_sectormb; i++) {
  			totalcharge += m_mbenergy[i];
  		}
  		//if (dataormc && central && totalcharge < centrality_bin[18]) { continue; }
  		if (dataormc && central) {
  			float max_charge = mc_centrality(max_cent);
  			float min_charge = mc_centrality(min_cent);
  			if (totalcharge > mc_centrality(min_cent) || totalcharge < mc_centrality(max_cent)) {
  				continue; 
  			}
  		}
  		if (!dataormc && central && (m_centbin > max_cent || m_centbin < min_cent)) { continue; }
  		h_mbd->Fill(totalcharge);

  		h_vz->Fill(m_vtx[2]);
  		if (reweighting && dataormc) { vz_weight = vertex_reweight[int(floor(m_vtx[2]*2)+100)]; }
  		h_vz_reweight->Fill(m_vtx[2],vz_weight);
  		totalweights += vz_weight;

  		if (dataormc) {
			for (int i = 0; i < m_g4; i++) {
	    		float theta = atan(m_g4_pt[i] / m_g4_pz[i]);
	    		float ET = m_g4_e[i] * abs(sin(theta));
	    		hetdeta->Fill(m_g4_eta[i], ET*vz_weight);
	    		hetdeta_emcalbin->Fill(m_g4_eta[i], ET*vz_weight);
	    		hetdeta_ihcalbin->Fill(m_g4_eta[i], ET*vz_weight);
	    		hetdeta_ohcalbin->Fill(m_g4_eta[i], ET*vz_weight);
	    		hetdeta_calobin->Fill(m_g4_eta[i], ET*vz_weight);
	    		if (fabs(m_g4_eta[i]) <= 1.1) {
	    			truthe += ET;
	    			hetdeta_zoom->Fill(m_g4_eta[i], ET*vz_weight);
	    		} 
			}
			h_event_truth_energy->Fill(truthe,vz_weight);
		}
		
		for (int i = 0; i < m_simtwrmult_cemc; i++) {
			//if (m_simtwr_cemc_ieta[i] < 8) { continue; }
			//std::tuple<int, int> hot_tower = std::make_tuple(m_simtwr_cemc_ieta[i], m_simtwr_cemc_iphi[i]);
		    //auto it = emcal_hot_dead_map.find(hot_tower);
		    //if (it != emcal_hot_dead_map.end()) { continue; }
			if (!dataormc && zs == 1 && m_simtwr_cemc_zs_adc[i] < zs_value) {
				float zs_e = 0; // = m_simtwr_cemc_zs_e[i]; 
				if (em_zs_calib[m_simtwr_cemc_ieta[i]][m_simtwr_cemc_iphi[i]] == 0) {
					h_em_zero_zscross->Fill(m_simtwr_cemc_ieta[i], m_simtwr_cemc_iphi[i]);
				} else {
					zs_e = m_simtwr_cemc_zs_e[i]/em_zs_calib[m_simtwr_cemc_ieta[i]][m_simtwr_cemc_iphi[i]];
				}
				h_2D_emcal_calib->Fill(m_simtwr_cemc_ieta[i], m_simtwr_cemc_iphi[i], zs_e*vz_weight);
				h_2D_emcal_calibT->Fill(m_simtwr_cemc_ieta[i], m_simtwr_cemc_iphi[i], zs_e*vz_weight/cosh(m_simtwr_cemc_eta[i]));
				emcale += zs_e/cosh(m_simtwr_cemc_eta[i]); 
				h_emcal->Fill(zs_e/cosh(m_simtwr_cemc_eta[i]), vz_weight);
				h_eT_emcal->Fill(m_simtwr_cemc_eta[i],zs_e*vz_weight/cosh(m_simtwr_cemc_eta[i]));
				h_eT_calo->Fill(m_simtwr_cemc_eta[i],zs_e*vz_weight/cosh(m_simtwr_cemc_eta[i]));
			} else if (!dataormc && zs == 2 && m_simtwr_cemc_zs_adc[i] < zs_value) {
				continue;
			} else {
				if (!dataormc && time && (m_simtwr_cemc_time[i] < -2 || m_simtwr_cemc_time[i] > 2)) { continue; }
				h_2D_emcal_calib->Fill(m_simtwr_cemc_ieta[i], m_simtwr_cemc_iphi[i], m_simtwr_cemc_e[i]*vz_weight);
				h_2D_emcal_calibT->Fill(m_simtwr_cemc_ieta[i], m_simtwr_cemc_iphi[i], m_simtwr_cemc_e[i]*vz_weight/cosh(m_simtwr_cemc_eta[i]));
				emcale += m_simtwr_cemc_e[i]/cosh(m_simtwr_cemc_eta[i]); 
				h_emcal->Fill(m_simtwr_cemc_e[i]/cosh(m_simtwr_cemc_eta[i]), vz_weight);
				h_eT_emcal->Fill(m_simtwr_cemc_eta[i],m_simtwr_cemc_e[i]*vz_weight/cosh(m_simtwr_cemc_eta[i]));
				h_eT_calo->Fill(m_simtwr_cemc_eta[i],m_simtwr_cemc_e[i]*vz_weight/cosh(m_simtwr_cemc_eta[i]));
			}
		}

		for (int i = 0; i < m_simtwrmult_ihcal; i++) {
			//std::tuple<int, int> hot_tower = std::make_tuple(m_simtwr_ihcal_ieta[i], m_simtwr_ihcal_iphi[i]);
		    //auto it = ihcal_hot_dead_map.find(hot_tower);
		    //if (it != ihcal_hot_dead_map.end()) { continue; }
		    if (!dataormc && zs == 1 && m_simtwr_ihcal_zs_adc[i] < zs_value) { 
		    	float zs_e = 0; // m_simtwr_ihcal_zs_e[i]; 
				if (ih_zs_calib[m_simtwr_ihcal_ieta[i]][m_simtwr_ihcal_iphi[i]] == 0) {
					h_ih_zero_zscross->Fill(m_simtwr_ihcal_ieta[i], m_simtwr_ihcal_iphi[i]);
				} else {
					zs_e = m_simtwr_ihcal_zs_e[i]/ih_zs_calib[m_simtwr_ihcal_ieta[i]][m_simtwr_ihcal_iphi[i]];
				}
				h_2D_ihcal_calib->Fill(m_simtwr_ihcal_ieta[i], m_simtwr_ihcal_iphi[i], zs_e*vz_weight);
				h_2D_ihcal_calibT->Fill(m_simtwr_ihcal_ieta[i], m_simtwr_ihcal_iphi[i], zs_e*vz_weight/cosh(m_simtwr_ihcal_eta[i]));
				ihcale += zs_e/cosh(m_simtwr_ihcal_eta[i]); 
				h_ihcal->Fill(zs_e/cosh(m_simtwr_ihcal_eta[i]), vz_weight);
				h_eT_ihcal->Fill(m_simtwr_ihcal_eta[i],zs_e*vz_weight/cosh(m_simtwr_ihcal_eta[i]));
				h_eT_calo->Fill(m_simtwr_ihcal_eta[i],zs_e*vz_weight/cosh(m_simtwr_ihcal_eta[i]));
		    } else if (!dataormc && zs == 2 && m_simtwr_ihcal_zs_adc[i] < zs_value) {
		    	continue;
			} else {
				h_2D_ihcal_calib->Fill(m_simtwr_ihcal_ieta[i], m_simtwr_ihcal_iphi[i], m_simtwr_ihcal_e[i]*vz_weight);
				h_2D_ihcal_calibT->Fill(m_simtwr_ihcal_ieta[i], m_simtwr_ihcal_iphi[i], m_simtwr_ihcal_e[i]*vz_weight/cosh(m_simtwr_ihcal_eta[i]));
				ihcale += m_simtwr_ihcal_e[i]/cosh(m_simtwr_ihcal_eta[i]); 
				h_ihcal->Fill(m_simtwr_ihcal_e[i]/cosh(m_simtwr_ihcal_eta[i]), vz_weight);
				h_eT_ihcal->Fill(m_simtwr_ihcal_eta[i],m_simtwr_ihcal_e[i]*vz_weight/cosh(m_simtwr_ihcal_eta[i]));
				h_eT_calo->Fill(m_simtwr_ihcal_eta[i],m_simtwr_ihcal_e[i]*vz_weight/cosh(m_simtwr_ihcal_eta[i]));
			}
		}

		for (int i = 0; i < m_simtwrmult_ohcal; i++) {
			//std::tuple<int, int> hot_tower = std::make_tuple(m_simtwr_ohcal_ieta[i], m_simtwr_ohcal_iphi[i]);
		    //auto it = ohcal_hot_dead_map.find(hot_tower);
		    //if (it != ohcal_hot_dead_map.end()) { continue; }
		    if (!dataormc && zs == 1 && m_simtwr_ohcal_zs_adc[i] < zs_value) {
		    	float zs_e = 0; // m_simtwr_ihcal_zs_e[i]; 
				if (oh_zs_calib[m_simtwr_ohcal_ieta[i]][m_simtwr_ohcal_iphi[i]] == 0) {
					h_oh_zero_zscross->Fill(m_simtwr_ohcal_ieta[i], m_simtwr_ohcal_iphi[i]);
				} else {
					zs_e = m_simtwr_ohcal_zs_e[i]/oh_zs_calib[m_simtwr_ohcal_ieta[i]][m_simtwr_ohcal_iphi[i]];
				}
		    	h_2D_ohcal_calib->Fill(m_simtwr_ohcal_ieta[i], m_simtwr_ohcal_iphi[i], zs_e*vz_weight);
				h_2D_ohcal_calibT->Fill(m_simtwr_ohcal_ieta[i], m_simtwr_ohcal_iphi[i], zs_e*vz_weight/cosh(m_simtwr_ohcal_eta[i]));
				ohcale += zs_e/cosh(m_simtwr_ohcal_eta[i]); 
				h_ohcal->Fill(zs_e/cosh(m_simtwr_ohcal_eta[i]), vz_weight);
				h_eT_ohcal->Fill(m_simtwr_ohcal_eta[i],zs_e*vz_weight/cosh(m_simtwr_ohcal_eta[i]));
				h_eT_calo->Fill(m_simtwr_ohcal_eta[i],zs_e*vz_weight/cosh(m_simtwr_ohcal_eta[i]));
		    } else if (!dataormc && zs == 2 && m_simtwr_ohcal_zs_adc[i] < zs_value) {
		    	continue;
			} else {
				h_2D_ohcal_calib->Fill(m_simtwr_ohcal_ieta[i], m_simtwr_ohcal_iphi[i], m_simtwr_ohcal_e[i]*vz_weight);
				h_2D_ohcal_calibT->Fill(m_simtwr_ohcal_ieta[i], m_simtwr_ohcal_iphi[i], m_simtwr_ohcal_e[i]*vz_weight/cosh(m_simtwr_ohcal_eta[i]));
				ohcale += m_simtwr_ohcal_e[i]/cosh(m_simtwr_ohcal_eta[i]); 
				h_ohcal->Fill(m_simtwr_ohcal_e[i]/cosh(m_simtwr_ohcal_eta[i]), vz_weight);
				h_eT_ohcal->Fill(m_simtwr_ohcal_eta[i],m_simtwr_ohcal_e[i]*vz_weight/cosh(m_simtwr_ohcal_eta[i]));
				h_eT_calo->Fill(m_simtwr_ohcal_eta[i],m_simtwr_ohcal_e[i]*vz_weight/cosh(m_simtwr_ohcal_eta[i]));
			}
		}
		
		totale = ihcale + ohcale;
		h_event_energy->Fill(totale + emcale, vz_weight);
		h_event_hcal_energy->Fill(totale, vz_weight);
		h_event_emcal_energy->Fill(emcale, vz_weight);
		h_event_ihcal_energy->Fill(ihcale, vz_weight);
		h_event_ohcal_energy->Fill(ohcale, vz_weight);
		
	}

	if (dataormc) {
		hetdeta->Scale(1.0/totalweights);
		hetdeta_zoom->Scale(1.0/totalweights);
		hetdeta->Scale(1.0/0.1);
		hetdeta_zoom->Scale(1.0/0.01);
		h_event_truth_energy->Scale(1.0/totalweights);
		hetdeta_ihcalbin->Scale(1.0/totalweights);
		hetdeta_ohcalbin->Scale(1.0/totalweights);
		hetdeta_emcalbin->Scale(1.0/totalweights);
		hetdeta_calobin->Scale(1.0/totalweights);
	}

	h_2D_emcal_calib->Scale(1.0/totalweights);
	h_2D_emcal_calibT->Scale(1.0/totalweights);
	h_emcal->Scale(1.0/totalweights);
	h_eT_emcal->Scale(1.0/totalweights);
	for (int i = 1; i <= emcal_num_bins; ++i) {
        double bin_width = emcal_bin_edges[i] - emcal_bin_edges[i - 1];
        h_eT_emcal->SetBinContent(i, h_eT_emcal->GetBinContent(i) / bin_width);
        h_eT_emcal->SetBinError(i, h_eT_emcal->GetBinError(i) / bin_width);
        if (dataormc) {
        	hetdeta_emcalbin->SetBinContent(i, hetdeta_emcalbin->GetBinContent(i)/ bin_width);
        	hetdeta_emcalbin->SetBinError(i, hetdeta_emcalbin->GetBinError(i)/ bin_width);
        }
    }

	h_2D_ihcal_calib->Scale(1.0/totalweights);
	h_2D_ihcal_calibT->Scale(1.0/totalweights);
	h_ihcal->Scale(1.0/totalweights);
	h_eT_ihcal->Scale(1.0/totalweights);
	for (int i = 1; i <= ihcal_num_bins; ++i) {
        double bin_width = ihcal_bin_edges[i] - ihcal_bin_edges[i - 1];
        h_eT_ihcal->SetBinContent(i, h_eT_ihcal->GetBinContent(i) / bin_width);
        h_eT_ihcal->SetBinError(i, h_eT_ihcal->GetBinError(i) / bin_width);
        if (dataormc) {
        	hetdeta_ihcalbin->SetBinContent(i, hetdeta_ihcalbin->GetBinContent(i)/ bin_width);
        	hetdeta_ihcalbin->SetBinError(i, hetdeta_ihcalbin->GetBinError(i)/ bin_width);
        }
    }

	h_2D_ohcal_calib->Scale(1.0/totalweights);
	h_2D_ohcal_calibT->Scale(1.0/totalweights);
	h_ohcal->Scale(1.0/totalweights);
	h_eT_ohcal->Scale(1.0/totalweights);
	for (int i = 1; i <= ohcal_num_bins; ++i) {
        double bin_width = ohcal_bin_edges[i] - ohcal_bin_edges[i - 1];
        h_eT_ohcal->SetBinContent(i, h_eT_ohcal->GetBinContent(i) / bin_width);
        h_eT_ohcal->SetBinError(i, h_eT_ohcal->GetBinError(i) / bin_width);
        if (dataormc) {
        	hetdeta_ohcalbin->SetBinContent(i, hetdeta_ohcalbin->GetBinContent(i) / bin_width);
        	hetdeta_ohcalbin->SetBinError(i, hetdeta_ohcalbin->GetBinError(i) / bin_width);
        }
    }

    h_eT_calo->Scale(1.0/totalweights);
	for (int i = 1; i <= calo_num_bins; ++i) {
        double bin_width = calo_bin_edges[i] - calo_bin_edges[i - 1];
        h_eT_calo->SetBinContent(i, h_eT_calo->GetBinContent(i) / bin_width);
        h_eT_calo->SetBinError(i, h_eT_calo->GetBinError(i) / bin_width);
        if (dataormc) {
        	hetdeta_calobin->SetBinContent(i, hetdeta_calobin->GetBinContent(i) / bin_width);
        	hetdeta_calobin->SetBinError(i, hetdeta_calobin->GetBinError(i) / bin_width);
        }
    }

	h_event_energy->Scale(1.0/totalweights);
	h_event_hcal_energy->Scale(1.0/totalweights);
	h_event_emcal_energy->Scale(1.0/totalweights);
	h_event_ihcal_energy->Scale(1.0/totalweights);
	h_event_ohcal_energy->Scale(1.0/totalweights);

	TH1F* h_emcal_correction;
	TH1F* h_ihcal_correction;
	TH1F* h_ohcal_correction;
	TH1F* h_calo_correction;

	if (dataormc) {
		h_emcal_correction = dynamic_cast<TH1F *>(h_eT_emcal->Clone("h_emcal_correction"));
		h_emcal_correction->Divide(hetdeta_emcalbin);
		h_ihcal_correction = dynamic_cast<TH1F *>(h_eT_ihcal->Clone("h_ihcal_correction"));
		h_ihcal_correction->Divide(hetdeta_ihcalbin);
		h_ohcal_correction = dynamic_cast<TH1F *>(h_eT_ohcal->Clone("h_ohcal_correction"));
		h_ohcal_correction->Divide(hetdeta_ohcalbin);
		h_calo_correction = dynamic_cast<TH1F *>(h_eT_calo->Clone("h_calo_correction"));
		h_calo_correction->Divide(hetdeta_calobin);
	}
	
	out->Write();
	out->Close();


}
