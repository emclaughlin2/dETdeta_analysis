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

const float hijing_cent[] = {3.0,1088.0,1941.0,3020.0,4478.0,6418.0,8983.0,12221.0,16172.0,21010.0,26857.0,33720.0,41673.0,50856.0,61295.0,73376.0,87385.0,103226.0,121633.0,143473.0,250000.0};
const float epos_cent[] = {27.0,1216.0,2122.0,3341.0,5071.0,7311.0,9977.0,13405.0,17700.0,22914.0,29130.0,36602.0,45471.0,55881.0,67429.0,81229.0,97175.0,115926.0,138395.0,166584.0,250000.0};
const float ampt_cent[] = {1.0,1206.0,2228.0,3562.0,5316.0,7688.0,10525.0,14400.0,19223.0,24427.0,30656.0,38090.0,47382.0,57843.0,69679.0,83789.0,100250.0,121187.0,142123.0,165105.0,250000.0};

std::set<std::tuple<int, int>> emcal_hot_dead_map_21615 = {{9, 32}, {9, 33}, {9, 34}, {9, 35}, {9, 36}, {9, 37}, {9, 38}, {9, 39}, {9, 40}, {9, 41}, {9, 42}, {9, 43}, {9, 44}, {9, 45}, {9, 46}, {9, 47}, {9, 136}, {9, 137}, {9, 138}, {9, 139}, {9, 140}, {9, 141}, {9, 142}, {9, 143}, {9, 200}, {10, 32}, {10, 33}, {10, 34}, {10, 35}, {10, 36}, {10, 37}, {10, 38}, {10, 39}, {10, 40}, {10, 41}, {10, 42}, {10, 43}, {10, 44}, {10, 45}, {10, 46}, {10, 47}, {10, 136}, {10, 137}, {10, 138}, {10, 139}, {10, 140}, {10, 141}, {10, 142}, {10, 143}, {11, 32}, {11, 33}, {11, 34}, {11, 35}, {11, 36}, {11, 37}, {11, 38}, {11, 39}, {11, 40}, {11, 41}, {11, 42}, {11, 43}, {11, 44}, {11, 45}, {11, 46}, {11, 47}, {11, 136}, {11, 137}, {11, 138}, {11, 139}, {11, 140}, {11, 141}, {11, 142}, {11, 143}, {12, 32}, {12, 33}, {12, 34}, {12, 35}, {12, 36}, {12, 37}, {12, 38}, {12, 39}, {12, 40}, {12, 41}, {12, 42}, {12, 43}, {12, 44}, {12, 45}, {12, 46}, {12, 47}, {12, 136}, {12, 137}, {12, 138}, {12, 139}, {12, 140}, {12, 141}, {12, 142}, {12, 143}, {13, 32}, {13, 33}, {13, 34}, {13, 35}, {13, 36}, {13, 37}, {13, 38}, {13, 39}, {13, 40}, {13, 41}, {13, 42}, {13, 43}, {13, 44}, {13, 45}, {13, 46}, {13, 47}, {13, 136}, {13, 137}, {13, 138}, {13, 139}, {13, 140}, {13, 141}, {13, 142}, {13, 143}, {13, 232}, {14, 32}, {14, 33}, {14, 34}, {14, 35}, {14, 36}, {14, 37}, {14, 38}, {14, 39}, {14, 40}, {14, 41}, {14, 42}, {14, 43}, {14, 44}, {14, 45}, {14, 46}, {14, 47}, {14, 136}, {14, 137}, {14, 138}, {14, 139}, {14, 140}, {14, 141}, {14, 142}, {14, 143}, {15, 32}, {15, 33}, {15, 34}, {15, 35}, {15, 36}, {15, 37}, {15, 38}, {15, 39}, {15, 40}, {15, 41}, {15, 42}, {15, 43}, {15, 44}, {15, 45}, {15, 46}, {15, 47}, {15, 136}, {15, 137}, {15, 138}, {15, 139}, {15, 140}, {15, 141}, {15, 142}, {15, 143}, {15, 215}, {16, 16}, {16, 17}, {16, 18}, {16, 19}, {16, 20}, {16, 21}, {16, 22}, {16, 23}, {16, 32}, {16, 33}, {16, 34}, {16, 35}, {16, 36}, {16, 37}, {16, 38}, {16, 39}, {16, 136}, {16, 137}, {16, 138}, {16, 139}, {16, 140}, {16, 141}, {16, 142}, {16, 143}, {16, 144}, {16, 145}, {16, 146}, {16, 147}, {16, 148}, {16, 149}, {16, 150}, {16, 151}, {17, 16}, {17, 17}, {17, 18}, {17, 19}, {17, 20}, {17, 21}, {17, 22}, {17, 23}, {17, 32}, {17, 33}, {17, 34}, {17, 35}, {17, 36}, {17, 37}, {17, 38}, {17, 39}, {17, 136}, {17, 137}, {17, 138}, {17, 139}, {17, 140}, {17, 141}, {17, 142}, {17, 143}, {17, 144}, {17, 145}, {17, 146}, 
{17, 147}, {17, 148}, {17, 149}, {17, 150}, {17, 151}, {18, 16}, {18, 17}, {18, 18}, {18, 19}, {18, 20}, {18, 21}, {18, 22}, {18, 23}, {18, 32}, {18, 33}, {18, 34}, {18, 35}, {18, 36}, {18, 37}, {18, 38}, {18, 39}, {18, 136}, {18, 137}, {18, 138}, {18, 139}, {18, 140}, {18, 141}, {18, 142}, {18, 143}, {18, 144}, {18, 145}, {18, 146}, {18, 147}, {18, 148}, {18, 149}, {18, 150}, {18, 151}, {19, 16}, {19, 17}, {19, 18}, {19, 19}, {19, 20}, {19, 21}, {19, 22}, {19, 23}, {19, 32}, {19, 33}, {19, 34}, {19, 35}, {19, 36}, {19, 37}, {19, 38}, {19, 39}, {19, 136}, {19, 137}, {19, 138}, {19, 139}, {19, 140}, {19, 141}, {19, 142}, {19, 143}, {19, 144}, {19, 145}, {19, 146}, {19, 147}, {19, 148}, {19, 149}, {19, 150}, {19, 151}, {20, 16}, {20, 17}, {20, 18}, {20, 19}, {20, 20}, {20, 21}, {20, 22}, {20, 23}, {20, 32}, {20, 33}, {20, 34}, {20, 35}, {20, 36}, {20, 37}, {20, 38}, {20, 39}, {20, 136}, {20, 137}, {20, 138}, {20, 139}, {20, 140}, {20, 141}, {20, 142}, {20, 143}, {20, 144}, {20, 145}, {20, 146}, {20, 147}, {20, 148}, {20, 149}, {20, 150}, {20, 151}, {21, 16}, {21, 17}, {21, 18}, {21, 19}, {21, 20}, {21, 21}, {21, 22}, {21, 23}, {21, 32}, {21, 33}, {21, 34}, {21, 35}, {21, 36}, {21, 37}, {21, 38}, {21, 39}, {21, 136}, {21, 137}, {21, 138}, {21, 139}, {21, 140}, {21, 141}, {21, 142}, {21, 143}, {21, 144}, {21, 145}, {21, 146}, {21, 147}, {21, 148}, {21, 149}, {21, 150}, {21, 151}, {22, 16}, {22, 17}, {22, 18}, {22, 19}, {22, 20}, {22, 21}, {22, 22}, {22, 23}, {22, 32}, {22, 33}, {22, 34}, {22, 35}, {22, 36}, {22, 37}, {22, 38}, {22, 39}, {22, 92}, {22, 119}, {22, 136}, {22, 137}, {22, 138}, {22, 139}, {22, 140}, {22, 141}, {22, 142}, {22, 143}, {22, 144}, {22, 145}, {22, 146}, {22, 147}, {22, 148}, {22, 149}, {22, 150}, {22, 151}, {22, 157}, {23, 16}, {23, 17}, {23, 18}, {23, 19}, {23, 20}, {23, 21}, {23, 22}, {23, 23}, {23, 32}, {23, 33}, {23, 34}, {23, 35}, {23, 36}, {23, 37}, {23, 38}, {23, 39}, {23, 136}, {23, 137}, {23, 138}, {23, 139}, {23, 140}, {23, 141}, {23, 142}, {23, 143}, {23, 144}, {23, 145}, {23, 146}, {23, 147}, {23, 148}, {23, 149}, {23, 150}, {23, 151}, {24, 32}, {24, 33}, {24, 34}, {24, 35}, {24, 36}, {24, 37}, {24, 38}, {24, 39}, {24, 78}, {24, 136}, {24, 137}, {24, 138}, {24, 139}, {24, 140}, {24, 141}, {24, 142}, {24, 143}, {24, 144}, {24, 145}, {24, 146}, {24, 147}, {24, 148}, {24, 149}, {24, 150}, {24, 151}, {24, 156}, {24, 208}, 
{24, 209}, {24, 210}, {24, 211}, {24, 212}, {24, 213}, {24, 214}, {24, 215}, {25, 32}, {25, 33}, {25, 34}, {25, 35}, {25, 36}, {25, 37}, {25, 38}, {25, 39}, {25, 136}, {25, 137}, {25, 138}, {25, 139}, {25, 140}, {25, 141}, {25, 142}, {25, 143}, {25, 144}, {25, 145}, {25, 146}, {25, 147}, {25, 148}, {25, 149}, {25, 150}, {25, 151}, {25, 208}, {25, 209}, {25, 210}, {25, 211}, {25, 212}, {25, 213}, {25, 214}, {25, 215}, {26, 32}, {26, 33}, {26, 34}, {26, 35}, {26, 36}, {26, 37}, {26, 38}, {26, 39}, {26, 136}, {26, 137}, {26, 138}, {26, 139}, {26, 140}, {26, 141}, {26, 142}, {26, 143}, {26, 144}, {26, 145}, {26, 146}, {26, 147}, {26, 148}, {26, 149}, {26, 150}, {26, 151}, {26, 208}, {26, 209}, {26, 210}, {26, 211}, {26, 212}, {26, 213}, {26, 214}, {26, 215}, {27, 32}, {27, 33}, {27, 34}, {27, 35}, {27, 36}, {27, 37}, {27, 38}, {27, 39}, {27, 136}, {27, 137}, {27, 138}, {27, 139}, {27, 140}, {27, 141}, {27, 142}, {27, 143}, {27, 144}, {27, 145}, {27, 146}, {27, 147}, {27, 148}, {27, 149}, {27, 150}, {27, 151}, {27, 208}, {27, 209}, {27, 210}, {27, 211}, {27, 212}, {27, 213}, {27, 214}, {27, 215}, {27, 223}, {28, 32}, {28, 33}, {28, 34}, {28, 35}, {28, 36}, {28, 37}, {28, 38}, {28, 39}, {28, 136}, {28, 137}, {28, 138}, {28, 139}, {28, 140}, {28, 141}, {28, 142}, {28, 143}, {28, 144}, {28, 145}, {28, 146}, {28, 147}, {28, 148}, {28, 149}, {28, 150}, {28, 151}, {28, 208}, {28, 209}, {28, 210}, {28, 211}, {28, 212}, {28, 213}, {28, 214}, {28, 215}, {29, 32}, {29, 33}, {29, 34}, {29, 35}, {29, 36}, {29, 37}, {29, 38}, {29, 39}, {29, 136}, {29, 137}, {29, 138}, {29, 139}, {29, 140}, {29, 141}, {29, 142}, {29, 143}, {29, 144}, {29, 145}, {29, 146}, {29, 147}, {29, 148}, {29, 149}, {29, 150}, {29, 151}, {29, 208}, {29, 209}, {29, 210}, {29, 211}, {29, 212}, {29, 213}, {29, 214}, {29, 215}, {30, 32}, {30, 33}, {30, 34}, {30, 35}, {30, 36}, {30, 37}, {30, 38}, {30, 39}, {30, 136}, {30, 137}, {30, 138}, {30, 139}, {30, 140}, {30, 141}, {30, 142}, {30, 143}, {30, 144}, {30, 145}, {30, 146}, {30, 147}, {30, 148}, {30, 149}, {30, 150}, {30, 151}, {30, 182}, {30, 208}, {30, 209}, {30, 210}, {30, 211}, {30, 212}, {30, 213}, {30, 214}, {30, 215}, {31, 32}, {31, 33}, {31, 34}, {31, 35}, {31, 36}, {31, 37}, {31, 38}, {31, 39}, {31, 104}, {31, 136}, {31, 137}, {31, 138}, {31, 139}, {31, 140}, {31, 141}, {31, 142}, {31, 143}, {31, 144}, {31, 145}, {31, 146}, {31, 147}, {31, 148}, {31, 149}, {31, 150}, {31, 151}, {31, 208}, {31, 209}, {31, 210}, {31, 211}, {31, 212}, {31, 213}, {31, 214}, {31, 215}, {32, 32}, {32, 33}, {32, 34}, {32, 35}, {32, 36}, {32, 37}, {32, 38}, {32, 39}, {32, 144}, {32, 145}, {32, 146}, {32, 147}, {32, 148}, {32, 149}, {32, 150}, {32, 151}, {33, 32}, {33, 33}, {33, 34}, {33, 35}, {33, 36}, {33, 37}, {33, 38}, {33, 39}, {33, 144}, {33, 145}, {33, 146}, {33, 147}, {33, 148}, {33, 149}, {33, 150}, {33, 151}, {34, 32}, {34, 33}, 
{34, 34}, {34, 35}, {34, 36}, {34, 37}, {34, 38}, {34, 39}, {34, 144}, {34, 145}, {34, 146}, {34, 147}, {34, 148}, {34, 149}, {34, 150}, {34, 151}, {35, 32}, {35, 33}, {35, 34}, {35, 35}, {35, 36}, {35, 37}, {35, 38}, {35, 39}, {35, 144}, {35, 145}, {35, 146}, {35, 147}, {35, 148}, {35, 149}, {35, 150}, {35, 151}, {36, 32}, {36, 33}, {36, 34}, {36, 35}, {36, 36}, {36, 37}, {36, 38}, {36, 39}, {36, 144}, {36, 145}, {36, 146}, {36, 147}, {36, 148}, {36, 149}, {36, 150}, {36, 151}, {37, 0}, {37, 32}, {37, 33}, {37, 34}, {37, 35}, {37, 36}, {37, 37}, {37, 38}, {37, 39}, {37, 144}, {37, 145}, {37, 146}, {37, 147}, {37, 148}, {37, 149}, {37, 150}, {37, 151}, {38, 32}, {38, 33}, {38, 34}, {38, 35}, {38, 36}, {38, 37}, {38, 38}, {38, 39}, {38, 144}, {38, 145}, {38, 146}, {38, 147}, {38, 148}, {38, 149}, {38, 150}, {38, 151}, {38, 219}, {39, 32}, {39, 33}, {39, 34}, {39, 35}, {39, 36}, {39, 37}, {39, 38}, {39, 39}, {39, 115}, {39, 144}, {39, 145}, {39, 146}, {39, 147}, {39, 148}, {39, 149}, {39, 150}, {39, 151}, {40, 32}, {40, 33}, {40, 34}, {40, 35}, {40, 36}, {40, 37}, {40, 38}, {40, 39}, {40, 91}, {40, 116}, {40, 144}, {40, 145}, {40, 146}, {40, 147}, {40, 148}, {40, 149}, {40, 150}, {40, 151}, {41, 32}, {41, 33}, {41, 34}, {41, 35}, {41, 36}, {41, 37}, {41, 38}, {41, 39}, {41, 144}, {41, 145}, {41, 146}, {41, 147}, {41, 148}, {41, 149}, {41, 150}, {41, 151}, {42, 32}, {42, 33}, {42, 34}, {42, 35}, {42, 36}, {42, 37}, {42, 38}, {42, 39}, {42, 144}, {42, 145}, {42, 146}, {42, 147}, {42, 148}, {42, 149}, {42, 150}, {42, 151}, {43, 32}, {43, 33}, {43, 34}, {43, 35}, {43, 36}, {43, 37}, {43, 38}, {43, 39}, {43, 144}, {43, 145}, {43, 146}, {43, 147}, {43, 148}, {43, 149}, {43, 150}, {43, 151}, {44, 32}, {44, 33}, {44, 34}, {44, 35}, {44, 36}, {44, 37}, {44, 38}, {44, 39}, {44, 144}, {44, 145}, {44, 146}, {44, 147}, {44, 148}, {44, 149}, {44, 150}, {44, 151}, {45, 32}, {45, 33}, {45, 34}, {45, 35}, {45, 36}, {45, 37}, {45, 38}, {45, 39}, {45, 144}, {45, 145}, {45, 146}, {45, 147}, {45, 148}, {45, 149}, {45, 150}, {45, 151}, {45, 220}, {46, 32}, {46, 33}, {46, 34}, {46, 35}, {46, 36}, {46, 37}, {46, 38}, {46, 39}, {46, 144}, {46, 145}, {46, 146}, {46, 147}, {46, 148}, {46, 149}, {46, 150}, {46, 151}, {47, 32}, {47, 33}, {47, 34}, {47, 35}, {47, 36}, {47, 37}, {47, 38}, {47, 39}, {47, 80}, {47, 96}, {47, 138}, {47, 144}, {47, 145}, {47, 146}, {47, 147}, {47, 148}, 
{47, 149}, {47, 150}, {47, 151}, {48, 31}, {48, 62}, {48, 79}, {48, 135}, {48, 139}, {48, 215}, {48, 231}, {48, 253}, {48, 255}, {49, 254}, {50, 243}, {52, 240}, {52, 246}, {55, 247}, {56, 160}, {57, 244}, {59, 241}, {63, 77}, {64, 88}, {64, 89}, {64, 90}, {64, 91}, {64, 92}, {64, 93}, {64, 94}, {64, 95}, {64, 164}, {65, 88}, {65, 89}, {65, 90}, {65, 91}, {65, 92}, {65, 93}, {65, 94}, {65, 95}, {65, 165}, {66, 88}, {66, 89}, {66, 90}, {66, 91}, {66, 92}, {66, 93}, {66, 94}, {66, 95}, {66, 246}, {67, 59}, {67, 88}, {67, 89}, {67, 90}, {67, 91}, {67, 92}, {67, 93}, {67, 94}, {67, 95}, {68, 88}, {68, 89}, {68, 90}, {68, 91}, {68, 92}, {68, 93}, {68, 94}, {68, 95}, {68, 139}, {68, 141}, {69, 88}, {69, 89}, {69, 90}, {69, 91}, {69, 92}, {69, 93}, {69, 94}, {69, 95}, {70, 88}, {70, 89}, {70, 90}, {70, 91}, {70, 92}, {70, 93}, {70, 94}, {70, 95}, {71, 88}, {71, 89}, {71, 90}, {71, 91}, {71, 92}, {71, 93}, {71, 94}, {71, 95}, {72, 188}, {72, 240}, {73, 189}, {78, 219}, {80, 149}, {80, 228}, {81, 150}, {85, 28}, {86, 27}, {86, 28}, {86, 173}, {88, 168}, {88, 169}, {88, 170}, {88, 171}, {88, 224}, {88, 225}, {88, 226}, {88, 227}, {89, 168}, {89, 169}, {89, 170}, {89, 171}, {89, 173}, {89, 224}, {89, 225}, {89, 226}, {89, 227}, {90, 111}, {90, 168}, {90, 169}, {90, 170}, {90, 171}, {90, 224}, {90, 225}, {90, 226}, {90, 227}, {91, 26}, {91, 151}, {91, 168}, {91, 169}, {91, 170}, {91, 171}, {91, 224}, {91, 225}, {91, 226}, {91, 227}, {91, 246}, {92, 0}, {92, 168}, {92, 169}, {92, 170}, {92, 171}, {92, 224}, {92, 225}, {92, 226}, {92, 227}, {93, 57}, {93, 96}, {93, 168}, {93, 169}, {93, 170}, {93, 171}, {93, 224}, {93, 225}, {93, 226}, {93, 227}, {94, 96}, {94, 106}, {94, 108}, {94, 157}, {94, 158}, {94, 159}, {94, 168}, {94, 169}, {94, 170}, {94, 171}, {94, 197}, {94, 219}, {94, 224}, {94, 225}, {94, 226}, {94, 227}, {95, 0}, {95, 1}, {95, 2}, {95, 3}, {95, 4}, {95, 5}, {95, 6}, {95, 7}, {95, 8}, {95, 9}, {95, 10}, {95, 11}, {95, 12}, {95, 13}, {95, 14}, {95, 15}, {95, 16}, {95, 17}, {95, 18}, {95, 19}, {95, 20}, {95, 21}, {95, 22}, {95, 23}, {95, 24}, {95, 25}, {95, 26}, {95, 27}, {95, 28}, {95, 29}, {95, 30}, {95, 31}, {95, 32}, {95, 33}, {95, 34}, {95, 35}, {95, 36}, {95, 37}, {95, 38}, {95, 39}, {95, 40}, {95, 41}, {95, 42}, {95, 43}, {95, 44}, {95, 45}, {95, 46}, {95, 47}, {95, 48}, {95, 49}, {95, 50}, {95, 51}, {95, 52}, {95, 53}, {95, 54}, {95, 55}, {95, 56}, 
{95, 57}, {95, 58}, {95, 59}, {95, 60}, {95, 61}, {95, 62}, {95, 63}, {95, 64}, {95, 65}, {95, 66}, {95, 67}, {95, 68}, {95, 69}, {95, 70}, {95, 71}, {95, 72}, {95, 73}, {95, 74}, {95, 75}, {95, 76}, {95, 77}, {95, 78}, {95, 79}, {95, 80}, {95, 81}, {95, 82}, {95, 83}, {95, 84}, {95, 85}, {95, 86}, {95, 87}, {95, 88}, {95, 89}, {95, 90}, {95, 91}, {95, 92}, {95, 93}, {95, 94}, {95, 95}, {95, 96}, {95, 97}, {95, 98}, {95, 99}, {95, 100}, {95, 101}, {95, 102}, {95, 103}, {95, 104}, {95, 105}, {95, 106}, {95, 107}, {95, 108}, {95, 109}, {95, 110}, {95, 111}, {95, 112}, {95, 113}, {95, 114}, {95, 115}, {95, 116}, {95, 117}, {95, 118}, {95, 119}, {95, 120}, {95, 121}, {95, 122}, {95, 123}, {95, 124}, {95, 125}, {95, 126}, {95, 127}, {95, 128}, {95, 129}, {95, 130}, {95, 131}, {95, 132}, {95, 133}, {95, 134}, {95, 135}, {95, 136}, {95, 137}, {95, 138}, {95, 139}, {95, 140}, {95, 141}, {95, 142}, {95, 143}, {95, 144}, {95, 145}, {95, 146}, {95, 147}, {95, 148}, {95, 149}, {95, 150}, {95, 151}, {95, 152}, {95, 153}, {95, 154}, {95, 155}, {95, 156}, {95, 157}, {95, 158}, {95, 159}, {95, 160}, {95, 161}, {95, 162}, {95, 163}, {95, 164}, {95, 165}, {95, 166}, {95, 167}, {95, 168}, {95, 169}, {95, 170}, {95, 171}, {95, 172}, {95, 173}, {95, 174}, {95, 175}, {95, 176}, {95, 177}, {95, 178}, {95, 179}, {95, 180}, {95, 181}, {95, 182}, {95, 183}, {95, 184}, {95, 185}, {95, 186}, {95, 187}, {95, 188}, {95, 189}, {95, 190}, {95, 191}, {95, 192}, {95, 193}, {95, 194}, {95, 195}, {95, 196}, {95, 197}, {95, 198}, {95, 199}, {95, 200}, {95, 201}, {95, 202}, {95, 203}, {95, 204}, {95, 205}, {95, 206}, {95, 207}, {95, 208}, {95, 209}, {95, 210}, {95, 211}, {95, 212}, {95, 213}, {95, 214}, {95, 215}, {95, 216}, {95, 217}, {95, 218}, {95, 219}, {95, 220}, {95, 221}, {95, 222}, {95, 223}, {95, 224}, {95, 225}, {95, 226}, {95, 227}, {95, 228}, {95, 229}, {95, 230}, {95, 231}, {95, 232}, {95, 233}, {95, 234}, {95, 235}, {95, 236}, {95, 237}, {95, 238}, {95, 239}, {95, 240}, {95, 241}, {95, 242}, {95, 243}, {95, 244}, {95, 245}, {95, 246}, {95, 247}, {95, 248}, {95, 249}, {95, 250}, {95, 251}, {95, 252}, {95, 253}, {95, 254}, {95, 255}};

const double vz_reweight[104] = {0.057359, 0.053709, 0.056632, 0.056701, 0.055505, 0.059661, 0.066884, 0.069031, 0.071042, 0.076489, 0.084428, 0.093634, 0.099436, 0.109863, 0.120297, 0.140174, 0.152318, 
0.175020, 0.187300, 0.208426, 0.247069, 0.287019, 0.338512, 0.369199, 0.435301, 0.511208, 0.565045, 0.667165, 0.748751, 0.913713, 1.091876, 1.170587, 1.422476, 1.633687, 1.919074, 2.204278, 2.553158, 
3.066737, 3.524783, 4.147713, 4.567779, 5.316948, 6.194884, 6.941905, 8.293592, 9.419249, 10.618740, 12.458636, 14.515464, 16.079992, 19.371290, 21.430447, 24.916803, 27.053852, 31.571571, 35.734322, 
40.051029, 43.771847, 47.863197, 50.340599, 55.608479, 56.607864, 66.653816, 67.427414, 69.564995, 71.192192, 81.833168, 78.077072, 83.215286, 82.865746, 85.609634, 81.214775, 82.856880, 87.721825, 
79.219582, 79.127502, 90.538330, 74.398285, 81.049438, 80.445770, 74.512970, 69.031090, 72.066505, 61.350029, 63.424370, 63.932232, 58.786026, 51.576374, 57.166645, 56.178078, 53.608730, 51.454712, 
46.199192, 46.557480, 41.882732, 48.681019, 43.275482, 41.458065, 41.662292, 29.244041, 35.175400, 36.495422, 37.713139, 36.780411};

const double ihcal_eta_bin_centers[24] = {-1.1419909849882701, -1.0565614038095177, -0.9712479248600292, -0.8859634212338496,
                            -0.8006991098291211, -0.7152269837198755, -0.6295610271831178, -0.5434190764827459,
                            -0.45669871862449246, -0.36929840324213237, -0.2810778018902582, -0.19168246028803004,
                            -0.10120418190943332, -0.009642418498081003, 0.0832040068717529, 0.17732534092707536,
                            0.27260507046668003, 0.36885741932987626, 0.4660686810181442, 0.5638669891904745,
                            0.6623151848090518, 0.7608530846929668, 0.8595314530240806, 0.958158998096158};

const double ohcal_eta_bin_centers[24] = {-1.1046054900559796, -1.0165613332375185, -0.9285735980639765, -0.8405913839757455, 
							-0.7525862447022397, -0.6644457661723757, -0.576116833503442, -0.4875010471303593, 
							-0.3985129668633035, -0.30908144991101183, -0.2190915218362084, -0.12847739861811402, 
							-0.03719644463345355, 0.054755588109315106, 0.14741290460362985, 0.24073652258128253, 
							0.33464084286045914, 0.4290447364594139, 0.5238869037463961, 0.6190367037577293, 
							0.714449171897191, 0.8099553357727513, 0.9054824008252557, 1.0009871856472095};

const int cemcSize = 24576;
const int ihcalSize = 1536;
const int ohcalSize = 1536;
const int mbdSize = 256;
const int g4Size = 100000;
const int vtxSize = 3;

void dETdeta_analysis(int dataormc = 0, int reweighting = 1, int central = 0) {	

	TDatabasePDG *_pdg = new TDatabasePDG();

	string filename = "dETdeta_analysis";
  	string dattag = (dataormc?"mc":"data");
  	string weighttag = (reweighting?"reweight":"noweight");
  	string centtag = (central?"0-10":"0-90");
  	filename += "_" + dattag + "_" + weighttag + "_" + centtag + "_epos.root";

	TFile *out = new TFile(filename.c_str(),"RECREATE");
	TH1F* h_vz = new TH1F("h_vz","",400, -100, 100);
	TH1F* h_vz_reweight = new TH1F("h_vz_reweight","",400, -100, 100);
	TH1F* h_mbd = new TH1F("h_mbd","",250000,0,250000);

	TH1F* h_event_truth_energy = new TH1F("h_event_truth_energy","",10000,0,10000);
	TH1F* hetdeta = new TH1F("hetdeta","",120,-6,6);
	TH1F* hetdeta_zoom = new TH1F("hetdeta_zoom","",220,-1.1,1.1);
	
	TH2F* h_2D_ihcal_calib = new TH2F("h_2D_ihcal_calib","",24,0.,24.,64,0.,64.);
  	TH2F* h_2D_ohcal_calib = new TH2F("h_2D_ohcal_calib","",24,0.,24.,64,0.,64.);
  	TH2F* h_2D_emcal_calib = new TH2F("h_2D_emcal_calib","",96,0.,96.,256,0.,256.);

  	TH2F* h_2D_ihcal_calibT = new TH2F("h_2D_ihcal_calibT","",24,0.,24.,64,0.,64.);
  	TH2F* h_2D_ohcal_calibT = new TH2F("h_2D_ohcal_calibT","",24,0.,24.,64,0.,64.);
  	TH2F* h_2D_emcal_calibT = new TH2F("h_2D_emcal_calibT","",96,0.,96.,256,0.,256.);
	
	TH1F* h_event_energy = new TH1F("h_event_energy","", 5000,0,5000);
	TH1F* h_event_hcal_energy = new TH1F("h_event_hcal_energy","", 5000,0,5000);
	TH1F* h_event_emcal_energy = new TH1F("h_event_emcal_energy","", 5000,0,5000);
	TH1F* h_event_ihcal_energy = new TH1F("h_event_ihcal_energy","", 1000,0,1000);
	TH1F* h_event_ohcal_energy = new TH1F("h_event_ohcal_energy","", 1000,0,1000);

	TH1F* h_emcal = new TH1F("h_emcal","",1000,0,10);
	TH1F* h_ihcal = new TH1F("h_ihcal","",1000,0,10);
	TH1F* h_ohcal = new TH1F("h_ohcal","",1000,0,10);

	TH1F* h_eT_emcal = new TH1F("h_eT_emcal","",220,-1.1,1.1);

	int ihcal_num_bins = 24;
    double ihcal_bin_edges[25];
    ihcal_bin_edges[0] = ihcal_eta_bin_centers[0] - 0.5 * (ihcal_eta_bin_centers[1] - ihcal_eta_bin_centers[0]);
    for (int i = 1; i < ihcal_num_bins; ++i) { ihcal_bin_edges[i] = (ihcal_eta_bin_centers[i] + ihcal_eta_bin_centers[i - 1]) / 2.0; }
    ihcal_bin_edges[ihcal_num_bins] = ihcal_eta_bin_centers[ihcal_num_bins - 1] + 0.5 * (ihcal_eta_bin_centers[ihcal_num_bins - 1] - ihcal_eta_bin_centers[ihcal_num_bins - 2]);
	TH1F* h_eT_ihcal = new TH1F("h_eT_ihcal","",ihcal_num_bins,ihcal_bin_edges);
	TH1F* h_eT_ihcal_noend = new TH1F("h_eT_ihcal_noend","",ihcal_num_bins,ihcal_bin_edges);
	TH1F* hetdeta_ihcalbin = new TH1F("hetdeta_ihcalbin","",ihcal_num_bins, ihcal_bin_edges);

	int ohcal_num_bins = 24;
    double ohcal_bin_edges[25];
    ohcal_bin_edges[0] = ohcal_eta_bin_centers[0] - 0.5 * (ohcal_eta_bin_centers[1] - ohcal_eta_bin_centers[0]);
    for (int i = 1; i < ohcal_num_bins; ++i) { ohcal_bin_edges[i] = (ohcal_eta_bin_centers[i] + ohcal_eta_bin_centers[i - 1]) / 2.0; }
    ohcal_bin_edges[ohcal_num_bins] = ohcal_eta_bin_centers[ohcal_num_bins - 1] + 0.5 * (ohcal_eta_bin_centers[ohcal_num_bins - 1] - ohcal_eta_bin_centers[ohcal_num_bins - 2]);
	TH1F* h_eT_ohcal = new TH1F("h_eT_ohcal","",ohcal_num_bins,ohcal_bin_edges);
	TH1F* h_eT_ohcal_noend = new TH1F("h_eT_ohcal_noend","",ohcal_num_bins,ohcal_bin_edges);
	TH1F* hetdeta_ohcalbin = new TH1F("hetdeta_ohcalbin","",ohcal_num_bins, ohcal_bin_edges);

	TH1F* h_eta_spread_ihcal[24];
	TH1F* h_eta_spread_ohcal[24];
	for (int i = 0; i < 24; i++) {
		h_eta_spread_ihcal[i] = new TH1F(TString::Format("h_eta_spread_ihcal_%d",i),"",400,-2,2);
		h_eta_spread_ohcal[i] = new TH1F(TString::Format("h_eta_spread_ohcal_%d",i),"",400,-2,2);
	}
    
    TChain chain("ttree");

    // location of hijing and data files
    //const char* inputDirectory = "/sphenix/user/egm2153/calib_study/detdeta/runsimana0/output/evt/";
	//TString wildcardPath = TString::Format("%sevents_20231122_nopileup_mc_cor*.root", inputDirectory);
    //if (!dataormc) { wildcardPath = TString::Format("%sevents_20231129_p004_data_cor*.root", inputDirectory); }
    //chain.Add(wildcardPath);

    // location of Joey's files
    //TString wildcardPath = "/sphenix/user/jocl/projects/sandbox/datatemp/merged_dEdeta_20231129_21615_mc_cor_555.root";
    //if (!dataormc) TString wildcardPath = "/sphenix/user/jocl/projects/sandbox/datatemp/merged_dEdeta_20231129_21615_data_cor_600.root";
    //chain.Add(wildcardPath);

    // location of EPOS files 
    const char* inputDirectory = "/sphenix/user/egm2153/calib_study/detdeta/eposrun/condor/";
    for (int i = 0; i < 900; i++) {
    	TString wildcardPath = TString::Format("%sOutDir%d/epos_sim_output.root", inputDirectory, i);
    	chain.Add(wildcardPath);
    }

    // testing files 
    //chain.Add("/sphenix/user/egm2153/calib_study/detdeta/runsimana0/output/evt/events_20231122_nopileup_mc_cor_0.root");
    //chain.Add("/sphenix/user/egm2153/calib_study/detdeta/runsimana0/output/evt/events_20231122_nopileup_mc_cor_1.root");
	//chain.Add("/sphenix/user/egm2153/calib_study/detdeta/runsimana0/output/evt/events_20231122_nopileup_mc_cor_2.root");
	//chain.Add("/sphenix/user/egm2153/calib_study/detdeta/runsimana0/output/evt/events_20231122_nopileup_mc_cor_3.root");
	//chain.Add("/sphenix/user/egm2153/calib_study/detdeta/runsimana0/output/evt/events_20231122_nopileup_mc_cor_4.root");

    int m_simtwrmult_cemc;
    float m_simtwr_cemc_e[cemcSize];
    int m_simtwr_cemc_ieta[cemcSize];
    int m_simtwr_cemc_iphi[cemcSize];
    float m_simtwr_cemc_eta[cemcSize];
    int m_simtwrmult_ihcal;
    float m_simtwr_ihcal_e[ihcalSize];
    int m_simtwr_ihcal_ieta[ihcalSize];
    int m_simtwr_ihcal_iphi[ihcalSize];
    float m_simtwr_ihcal_eta[ihcalSize];
    int m_simtwrmult_ohcal;
    float m_simtwr_ohcal_e[ohcalSize];
    int m_simtwr_ohcal_ieta[ohcalSize];
    int m_simtwr_ohcal_iphi[ohcalSize];
    float m_simtwr_ohcal_eta[ohcalSize];
    int m_sectormb;
    float m_mbenergy[mbdSize];
    int m_g4;
    float m_g4_e[g4Size];
    float m_g4_eta[g4Size];
    float m_g4_pt[g4Size];
    float m_g4_pz[g4Size];
    float m_vtx[vtxSize];
    float mbd_vtx[vtxSize];

     // Set branch addresses
    chain.SetBranchAddress("sectorem", &m_simtwrmult_cemc);
    chain.SetBranchAddress("emcalen", m_simtwr_cemc_e);
    chain.SetBranchAddress("emcaletabin", m_simtwr_cemc_ieta);
    chain.SetBranchAddress("emcalphibin", m_simtwr_cemc_iphi);
    chain.SetBranchAddress("emetacor", m_simtwr_cemc_eta);

    chain.SetBranchAddress("sectorih", &m_simtwrmult_ihcal);
    chain.SetBranchAddress("ihcalen", m_simtwr_ihcal_e);
    chain.SetBranchAddress("ihcaletabin", m_simtwr_ihcal_ieta);
    chain.SetBranchAddress("ihcalphibin", m_simtwr_ihcal_iphi);
    chain.SetBranchAddress("ihetacor", m_simtwr_ihcal_eta);

    chain.SetBranchAddress("sectoroh", &m_simtwrmult_ohcal);
    chain.SetBranchAddress("ohcalen", m_simtwr_ohcal_e);
    chain.SetBranchAddress("ohcaletabin", m_simtwr_ohcal_ieta);
    chain.SetBranchAddress("ohcalphibin", m_simtwr_ohcal_iphi);
    chain.SetBranchAddress("ohetacor", m_simtwr_ohcal_eta);

    chain.SetBranchAddress("sectormb", &m_sectormb);
    chain.SetBranchAddress("mbenrgy", &m_mbenergy);

    if (dataormc) {
    	chain.SetBranchAddress("truthpar_n", &m_g4);
	    chain.SetBranchAddress("truthpar_e", m_g4_e);
	    chain.SetBranchAddress("truthpar_eta", m_g4_eta);
	    chain.SetBranchAddress("truthpar_pt", m_g4_pt);
	    chain.SetBranchAddress("truthpar_pz", m_g4_pz);
    }

    chain.SetBranchAddress("track_vtx", m_vtx);
    chain.SetBranchAddress("mbd_vtx", mbd_vtx);

	int eventnumber = 0;
	float totalweights = 0.0;
	float delta_eta = 0.09167;
	float delta_em_eta = 0.022918;

	float data_centrality_binning[10] = {0.0,27.0,66.0,131.0,228.0,368.0,573.0,860.0,1247.0,250000.0};
	//float mc_centrality_binning[10] = {1915.0,4443.0,8920.0,16156.0,26795.0,41511.0,61118.0,86986.0,121515.0,250000.0}; // hijing
	float mc_centrality_binning[10] = {2254.0,5152.0,10507.0,17717.0,28138.0,43228.0,65105.0,93447.0,132042.0,250000.0}; // EPOS

    std::vector<float> centrality_bin;
    if (dataormc) centrality_bin.assign(mc_centrality_binning, mc_centrality_binning+10);
    if (!dataormc) centrality_bin.assign(data_centrality_binning, data_centrality_binning+10);

    assert(centrality_bin.size() == 10);

    Long64_t nEntries = chain.GetEntries();
    for (Long64_t entry = 0; entry < nEntries; ++entry) {
    //for (Long64_t entry = 0; entry < 100; ++entry) {
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
  		// require that simulation could reconstruct a vertex for the event
  		if (m_vtx[2] <= -6.0 || m_vtx[2] >= 46.0) { continue; }

  		float totalcharge = 0.0;
  		for (int i = 0; i < m_sectormb; i++) {
  			totalcharge += m_mbenergy[i];
  		}
  		if (central && totalcharge < centrality_bin[8]) { continue; }
  		h_mbd->Fill(totalcharge);

  		h_vz->Fill(m_vtx[2]);
  		if (reweighting) { vz_weight = vz_reweight[int(floor(m_vtx[2]*2)+12)]; }
  		h_vz_reweight->Fill(m_vtx[2],vz_weight);
  		totalweights += vz_weight;

  		if (dataormc) {
			for (int i = 0; i < m_g4; i++) {
	    		float theta = atan(m_g4_pt[i] / m_g4_pz[i]);
	    		float ET = m_g4_e[i] * abs(sin(theta));
	    		hetdeta->Fill(m_g4_eta[i], ET*vz_weight);
	    		hetdeta_ihcalbin->Fill(m_g4_eta[i], ET*vz_weight);
	    		hetdeta_ohcalbin->Fill(m_g4_eta[i], ET*vz_weight);
	    		if (fabs(m_g4_eta[i]) < 1.1) {
	    			truthe += ET;
	    			hetdeta_zoom->Fill(m_g4_eta[i], ET*vz_weight);
	    		} 
			}
			h_event_truth_energy->Fill(truthe,vz_weight);
		}
		
		for (int i = 0; i < m_simtwrmult_cemc; i++) {
			if (m_simtwr_cemc_ieta[i] < 9) { continue; }
			if (m_simtwr_cemc_ieta[i] >= 48 && m_simtwr_cemc_ieta[i] <= 95 && m_simtwr_cemc_iphi[i] >= 240 && m_simtwr_cemc_iphi[i] <= 247) { continue; }
			std::tuple<int, int> hot_tower = std::make_tuple(m_simtwr_cemc_ieta[i], m_simtwr_cemc_iphi[i]);
		    auto it = emcal_hot_dead_map_21615.find(hot_tower);
		    if (it != emcal_hot_dead_map_21615.end()) { continue; }
			h_2D_emcal_calib->Fill(m_simtwr_cemc_ieta[i], m_simtwr_cemc_iphi[i], m_simtwr_cemc_e[i]*vz_weight);
			h_2D_emcal_calibT->Fill(m_simtwr_cemc_ieta[i], m_simtwr_cemc_iphi[i], m_simtwr_cemc_e[i]*vz_weight/cosh(m_simtwr_cemc_eta[i]));
			emcale += m_simtwr_cemc_e[i]/cosh(m_simtwr_cemc_eta[i]); 
			h_emcal->Fill(m_simtwr_cemc_e[i]/cosh(m_simtwr_cemc_eta[i]), vz_weight);
			h_eT_emcal->Fill(m_simtwr_cemc_eta[i],m_simtwr_cemc_e[i]*vz_weight/cosh(m_simtwr_cemc_eta[i]));
		}


		for (int i = 0; i < m_simtwrmult_ihcal; i++) {
			if (m_simtwr_ihcal_ieta[i] >= 0 && m_simtwr_ihcal_ieta[i] <= 7 && m_simtwr_ihcal_iphi[i] >= 48 && m_simtwr_ihcal_iphi[i] <= 49) { continue; }
			if (m_simtwr_ihcal_ieta[i] == 21 && m_simtwr_ihcal_iphi[i] == 49) { continue; }
			if (m_simtwr_ihcal_ieta[i] == 8 && m_simtwr_ihcal_iphi[i] == 32) { continue; }
			if (m_simtwr_ihcal_ieta[i] == 7 && m_simtwr_ihcal_iphi[i] == 51) { continue; }
			h_2D_ihcal_calib->Fill(m_simtwr_ihcal_ieta[i], m_simtwr_ihcal_iphi[i], m_simtwr_ihcal_e[i]*vz_weight);
			h_2D_ihcal_calibT->Fill(m_simtwr_ihcal_ieta[i], m_simtwr_ihcal_iphi[i], m_simtwr_ihcal_e[i]*vz_weight/cosh(m_simtwr_ihcal_eta[i]));
			ihcale += m_simtwr_ihcal_e[i]/cosh(m_simtwr_ihcal_eta[i]); 
			h_ihcal->Fill(m_simtwr_ihcal_e[i]/cosh(m_simtwr_ihcal_eta[i]), vz_weight);
			h_eT_ihcal->Fill(m_simtwr_ihcal_eta[i],m_simtwr_ihcal_e[i]*vz_weight/cosh(m_simtwr_ihcal_eta[i]));
			if (m_simtwr_ihcal_ieta[i] >= 2 && m_simtwr_ihcal_ieta[i] <= 21) h_eT_ihcal_noend->Fill(m_simtwr_ihcal_eta[i],m_simtwr_ihcal_e[i]*vz_weight/cosh(m_simtwr_ihcal_eta[i]));
			h_eta_spread_ihcal[m_simtwr_ihcal_ieta[i]]->Fill(m_simtwr_ihcal_eta[i], vz_weight);
		}


		for (int i = 0; i < m_simtwrmult_ohcal; i++) {
		    if (m_simtwr_ohcal_ieta[i] == 7 && m_simtwr_ohcal_iphi[i] == 31) { continue; }
			h_2D_ohcal_calib->Fill(m_simtwr_ohcal_ieta[i], m_simtwr_ohcal_iphi[i], m_simtwr_ohcal_e[i]*vz_weight);
			h_2D_ohcal_calibT->Fill(m_simtwr_ohcal_ieta[i], m_simtwr_ohcal_iphi[i], m_simtwr_ohcal_e[i]*vz_weight/cosh(m_simtwr_ohcal_eta[i]));
			ohcale += m_simtwr_ohcal_e[i]/cosh(m_simtwr_ohcal_eta[i]); 
			h_ohcal->Fill(m_simtwr_ohcal_e[i]/cosh(m_simtwr_ohcal_eta[i]), vz_weight);
			h_eT_ohcal->Fill(m_simtwr_ohcal_eta[i],m_simtwr_ohcal_e[i]*vz_weight/cosh(m_simtwr_ohcal_eta[i]));
			if (m_simtwr_ohcal_ieta[i] >= 2 && m_simtwr_ohcal_ieta[i] <= 21) h_eT_ohcal_noend->Fill(m_simtwr_ohcal_eta[i],m_simtwr_ohcal_e[i]*vz_weight/cosh(m_simtwr_ohcal_eta[i]));
			h_eta_spread_ohcal[m_simtwr_ohcal_ieta[i]]->Fill(m_simtwr_ohcal_eta[i], vz_weight);
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
	}

	h_2D_emcal_calib->Scale(1.0/totalweights);
	h_2D_emcal_calibT->Scale(1.0/totalweights);
	h_emcal->Scale(1.0/totalweights);
	h_eT_emcal->Scale(1.0/totalweights);
	h_eT_emcal->Scale(1.0/0.01);

	h_2D_ihcal_calib->Scale(1.0/totalweights);
	h_2D_ihcal_calibT->Scale(1.0/totalweights);
	h_ihcal->Scale(1.0/totalweights);
	h_eT_ihcal->Scale(1.0/totalweights);
	h_eT_ihcal_noend->Scale(1.0/totalweights);
	for (int i = 1; i <= ihcal_num_bins; ++i) {
        double bin_width = ihcal_bin_edges[i] - ihcal_bin_edges[i - 1];
        h_eT_ihcal->SetBinContent(i, h_eT_ihcal->GetBinContent(i) / bin_width);
        h_eT_ihcal_noend->SetBinContent(i, h_eT_ihcal_noend->GetBinContent(i) / bin_width);
        if (dataormc) hetdeta_ihcalbin->SetBinContent(i, hetdeta_ihcalbin->GetBinContent(i)/ bin_width);
    }

	h_2D_ohcal_calib->Scale(1.0/totalweights);
	h_2D_ohcal_calibT->Scale(1.0/totalweights);
	h_ohcal->Scale(1.0/totalweights);
	h_eT_ohcal->Scale(1.0/totalweights);
	h_eT_ohcal_noend->Scale(1.0/totalweights);
	for (int i = 1; i <= ohcal_num_bins; ++i) {
        double bin_width = ohcal_bin_edges[i] - ohcal_bin_edges[i - 1];
        h_eT_ohcal->SetBinContent(i, h_eT_ohcal->GetBinContent(i) / bin_width);
        h_eT_ohcal_noend->SetBinContent(i, h_eT_ohcal_noend->GetBinContent(i) / bin_width);
        if (dataormc) hetdeta_ohcalbin->SetBinContent(i, hetdeta_ohcalbin->GetBinContent(i) / bin_width);
    }

	for (int i = 0; i < 24; i++) {
		h_eta_spread_ihcal[i]->Scale(1.0/totalweights);
		h_eta_spread_ohcal[i]->Scale(1.0/totalweights);
		h_eta_spread_ihcal[i]->Scale(1.0/0.01);
		h_eta_spread_ohcal[i]->Scale(1.0/0.01);
	}

	h_event_energy->Scale(1.0/totalweights);
	h_event_hcal_energy->Scale(1.0/totalweights);
	h_event_emcal_energy->Scale(1.0/totalweights);
	h_event_ihcal_energy->Scale(1.0/totalweights);
	h_event_ohcal_energy->Scale(1.0/totalweights);
	
	out->Write();
	out->Close();


}
