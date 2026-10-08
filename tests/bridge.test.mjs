// Tests of the external-solver bridge: node tests/bridge.test.mjs   (non-zero exit on failure)
// Hand-off table against the catalogue, every generator on the default case, mesh topology and volume, the zip writer,
// the importers on output written by the real solvers (excerpts below), and the property-table round trip.
// Optional: BRIDGE_PYTHON=/path/to/python (with CoolProp and teqp installed) also executes the generated property scripts;
// HANDOFF_REQUESTS=/path/to/folder checks request files ([{ item, why, solver }] per <suite>.json) against the table.
import fs from 'node:fs';
import os from 'node:os';
import path from 'node:path';
import { execFileSync } from 'node:child_process';
import * as B from '../js/core/bridge.js';
import { CATALOG } from '../js/data/catalog.js';
import { buildTable, lookup, fluidModel, DEFAULT_FLUID } from '../js/core/thermo.js';

let fails = 0, count = 0;
const ok = (cond, msg) => { count++; if (!cond) { fails++; console.log('   ✗ ' + msg); } };
const near = (a, b, rel, msg) => ok(Number.isFinite(a) && Math.abs(a - b) <= rel * Math.max(Math.abs(b), 1e-300), `${msg}: got ${a}, expected ${b} within ${rel * 100} %`);
const section = (name) => console.log('• ' + name);
const SUITE_IDS = ['pvt', 'net', 'flow', 'solids', 'ops', 'integ', 'econ'];
const itemsOf = (id) => CATALOG[SUITE_IDS.indexOf(id) + 1].groups.flatMap((g) => g.items);

// Excerpts of files written by the real solvers when the generated cases were run during development.
const SAMPLE = {
  // interFoam, k-omega SST: the generated pipe-section case run in OpenFOAM v2412 (first rows of five function-object files, as collect.sh gathers them)
  pipe: `### file: postProcessing/forces/0/force.dat
# Force         
# CofR          : (0.00000000e+00 0.00000000e+00 0.00000000e+00)
#
# Time          	total_x total_y total_z	pressure_x pressure_y pressure_z	viscous_x viscous_y viscous_z
0.00125          5.19941890e-01 5.55860161e-04 -4.72674034e+02 -7.65822363e-15 5.58080028e-04 -4.72768277e+02 5.19941890e-01 -2.21986637e-06 9.42433941e-02
0.0025           5.86085957e-01 8.15033310e-04 -4.17562538e+02 -7.23865704e-15 8.16929297e-04 -4.17666494e+02 5.86085957e-01 -1.89598747e-06 1.03956168e-01
0.00375          6.24404859e-01 3.57462876e-04 -3.97972594e+02 -6.64393960e-15 3.59364521e-04 -3.98079335e+02 6.24404859e-01 -1.90164486e-06 1.06741068e-01
0.005            6.53651014e-01 1.12462550e-03 -3.87892133e+02 -6.27338102e-15 1.12695931e-03 -3.88000261e+02 6.53651014e-01 -2.33381758e-06 1.08128153e-01
0.00625          6.79375401e-01 3.18243879e-04 -3.81680589e+02 -6.00524825e-15 3.19802557e-04 -3.81789682e+02 6.79375401e-01 -1.55867870e-06 1.09092894e-01
0.0075           7.03656830e-01 -1.12070258e-03 -3.77469239e+02 -5.78221920e-15 -1.12177415e-03 -3.77579211e+02 7.03656830e-01 1.07157032e-06 1.09971941e-01
0.00875          7.27730818e-01 -1.12366901e-03 -3.74358168e+02 -5.61290633e-15 -1.12756988e-03 -3.74469147e+02 7.27730818e-01 3.90087105e-06 1.10978749e-01
0.01             7.52955772e-01 -2.96660626e-03 -3.72098721e+02 -5.48812734e-15 -2.97172344e-03 -3.72210806e+02 7.52955772e-01 5.11718071e-06 1.12085520e-01
0.01125          7.78398000e-01 -2.57843429e-03 -3.70588037e+02 -5.39547312e-15 -2.58282421e-03 -3.70701428e+02 7.78398000e-01 4.38992725e-06 1.13391662e-01
0.0125           8.06189715e-01 -2.86344932e-03 -3.69417635e+02 -5.34305738e-15 -2.87045384e-03 -3.69532386e+02 8.06189715e-01 7.00451586e-06 1.14751304e-01
0.01375          8.38008094e-01 -2.89551444e-03 -3.68457415e+02 -5.32851374e-15 -2.90640145e-03 -3.68573572e+02 8.38008094e-01 1.08870065e-05 1.16156881e-01
0.015            8.74737897e-01 -3.71719414e-03 -3.67776082e+02 -5.34188687e-15 -3.73157878e-03 -3.67893758e+02 8.74737897e-01 1.43846374e-05 1.17675602e-01

### file: postProcessing/holdupVolume/0/volFieldValue.dat
# Region        : all region0
# Cells         : 13920
# Volume        : 7.52526948e-02
# Time          	volAverage(alpha.liquid)
0.00125         	5.47289657e-01
0.0025          	5.48013711e-01
0.00375         	5.48736799e-01
0.005           	5.49458884e-01
0.00625         	5.50179930e-01
0.0075          	5.50899901e-01
0.00875         	5.51618761e-01
0.01            	5.52336476e-01
0.01125         	5.53053014e-01
0.0125          	5.53768348e-01
0.01375         	5.54482452e-01
0.015           	5.55195303e-01

### file: postProcessing/probes/0/p
# Probe 0 (0.381 0.01778 0.01778)
# Probe 1 (0.762 0.01778 0.01778)
# Probe 2 (1.143 0.01778 0.01778)
# Time          0               1               2              
0.00125         117.19735       79.374604       32.795167      
0.0025          183.8986        124.4431        52.811309      
0.00375         243.45102       165.48757       71.216421      
0.005           298.00014       203.32915       88.208421      
0.00625         349.12153       238.67781       104.10024      
0.0075          398.03374       271.61613       118.9178       
0.00875         444.28531       302.03756       132.59273      
0.01            487.95056       329.86569       145.06295      
0.01125         527.92156       355.11038       156.30416      
0.0125          564.5005        378.00065       166.42523      
0.01375         597.80153       398.82387       175.57361      
0.015           628.21804       417.80566       183.84078      

### file: postProcessing/section1/0/surfaceFieldValue.dat
# Region type : sampledSurface section1
# Faces         : 464
# Area          : 4.93784086e-02
# Scale factor  : 1.00000000e+00
# Time          	areaAverage(alpha.liquid)	areaAverage(p)
0.00125         	5.46737969e-01	3.47541190e+02
0.0025          	5.46784728e-01	4.13016898e+02
0.00375         	5.46842743e-01	4.73763662e+02
0.005           	5.46912302e-01	5.31371134e+02
0.00625         	5.46994016e-01	5.85899090e+02
0.0075          	5.47088254e-01	6.36927156e+02
0.00875         	5.47194989e-01	6.83842405e+02
0.01            	5.47313890e-01	7.26537085e+02
0.01125         	5.47444505e-01	7.65106555e+02
0.0125          	5.47586324e-01	7.99724031e+02
0.01375         	5.47738958e-01	8.30962417e+02
0.015           	5.47902040e-01	8.59255158e+02

### file: postProcessing/section3/0/surfaceFieldValue.dat
# Region type : sampledSurface section3
# Faces         : 464
# Area          : 4.93784086e-02
# Scale factor  : 1.00000000e+00
# Time          	areaAverage(alpha.liquid)	areaAverage(p)
0.00125         	5.46709723e-01	2.32360271e+02
0.0025          	5.46709675e-01	2.52622320e+02
0.00375         	5.46709592e-01	2.71193142e+02
0.005           	5.46709474e-01	2.88419901e+02
0.00625         	5.46709324e-01	3.04584586e+02
0.0075          	5.46709142e-01	3.19710841e+02
0.00875         	5.46708928e-01	3.33721237e+02
0.01            	5.46708685e-01	3.46563502e+02
0.01125         	5.46708410e-01	3.58228986e+02
0.0125          	5.46708106e-01	3.68803605e+02
0.01375         	5.46707772e-01	3.78417432e+02
0.015           	5.46707409e-01	3.87157987e+02

`,
  // compressibleInterFoam: the generated 2-D pipeline-riser case run for 22 s in OpenFOAM v2412 (riser-base and riser-top sections, every 25th row)
  riser: `### file: postProcessing/section2/0/surfaceFieldValue.dat
# Region type : sampledSurface section2
# Faces         : 10
# Area          : 6.45160008e-02
# Scale factor  : 1.00000000e+00
# Time          	areaAverage(alpha.liquid)	areaAverage(p)
0.0031113797    	2.00150754e-41	8.56061660e+06
0.36329334      	1.07375328e-04	8.56234036e+06
0.72066189      	6.25318200e-01	8.55897196e+06
1.0798016       	7.73682057e-01	8.56136390e+06
1.3989534       	7.82035254e-01	8.56178068e+06
1.7042669       	7.90501703e-01	8.56352040e+06
2.0054099       	7.99685105e-01	8.56558282e+06
2.3105285       	8.05123887e-01	8.56693062e+06
2.6222961       	8.09863723e-01	8.56924431e+06
2.9409577       	8.11987937e-01	8.57090163e+06
3.3149079       	8.58908878e-01	8.57342235e+06
3.8034098       	8.99875457e-01	8.57941121e+06
4.3920279       	9.45694757e-01	8.58510619e+06
5.0183699       	8.86400144e-01	8.59228472e+06
5.5596276       	8.62567376e-01	8.59606629e+06
5.9155984       	8.66348978e-01	8.59759326e+06
6.2102353       	8.66077441e-01	8.59727752e+06
6.4926321       	8.64070478e-01	8.59507391e+06
6.7476079       	8.65084932e-01	8.59749188e+06
7.0082041       	8.64501141e-01	8.60224499e+06
7.2625014       	8.62104396e-01	8.59848791e+06
7.5019608       	8.62244097e-01	8.59865129e+06
7.7362922       	8.61796503e-01	8.59965721e+06
7.9856821       	8.60534228e-01	8.60080843e+06
8.2618065       	8.58711313e-01	8.59990157e+06
8.5415587       	8.58772541e-01	8.60259974e+06
8.7844041       	8.57618543e-01	8.60090165e+06
9.0300112       	8.57492937e-01	8.60467602e+06
9.2673174       	8.56029936e-01	8.60252065e+06
9.5113947       	8.55093501e-01	8.59945703e+06
9.7491943       	8.54449264e-01	8.59732950e+06
9.9847677       	8.55158133e-01	8.60132144e+06
10.228487       	8.54040573e-01	8.60201786e+06
10.473112       	8.53720038e-01	8.60442864e+06
10.715338       	8.52290535e-01	8.60193364e+06
10.962178       	8.51907444e-01	8.60086438e+06
11.203224       	8.51655136e-01	8.60122217e+06
11.45013        	8.50445019e-01	8.59834238e+06
11.69024        	8.51130971e-01	8.60014928e+06
11.929091       	8.50881556e-01	8.60148608e+06
12.170959       	8.50344870e-01	8.60128807e+06
12.418552       	8.49783654e-01	8.60143478e+06
12.659236       	8.50013079e-01	8.60156941e+06
12.898039       	8.49712015e-01	8.60294027e+06
13.141231       	8.49385288e-01	8.60296157e+06
13.383018       	8.49337031e-01	8.60228367e+06
13.624803       	8.48800321e-01	8.60032658e+06
13.870181       	8.48242238e-01	8.59786055e+06
14.117744       	8.48706161e-01	8.59978452e+06
14.35358        	8.49007295e-01	8.60172252e+06
14.59921        	8.48368339e-01	8.60074130e+06
14.837704       	8.48111383e-01	8.59847924e+06
15.080992       	8.48532955e-01	8.60058035e+06
15.326917       	8.48641984e-01	8.60394229e+06
15.566603       	8.48740312e-01	8.60452081e+06
15.803277       	8.47891683e-01	8.60134640e+06
16.048048       	8.48096238e-01	8.60122603e+06
16.292131       	8.48387764e-01	8.60237852e+06
16.532814       	8.47834056e-01	8.60106748e+06
16.776131       	8.47907260e-01	8.60036037e+06
17.020316       	8.47642315e-01	8.59849153e+06
17.255744       	8.48185349e-01	8.60077543e+06
17.496635       	8.47967330e-01	8.60104400e+06
17.734876       	8.47864029e-01	8.60017810e+06
17.979185       	8.47853232e-01	8.59993498e+06
18.222654       	8.47999588e-01	8.60163986e+06
18.470661       	8.48224205e-01	8.60316982e+06
18.716408       	8.48116224e-01	8.60358302e+06
18.965916       	8.47121140e-01	8.59939156e+06
19.21266        	8.48116989e-01	8.60168445e+06
19.453065       	8.47938972e-01	8.60125893e+06
19.692758       	8.47975446e-01	8.60156206e+06
19.941695       	8.47177379e-01	8.59887384e+06
20.188445       	8.47929003e-01	8.60097948e+06
20.433951       	8.48048838e-01	8.60236511e+06
20.679175       	8.47695464e-01	8.60138634e+06
20.929055       	8.48096121e-01	8.60235606e+06
21.176505       	8.47295706e-01	8.59929025e+06
21.415867       	8.48065540e-01	8.60022006e+06
21.663796       	8.47665732e-01	8.59972613e+06
21.909817       	8.47921260e-01	8.59963388e+06
22.156909       	8.48486265e-01	8.60261312e+06
22.403709       	8.47615242e-01	8.60007681e+06

### file: postProcessing/section3/0/surfaceFieldValue.dat
# Region type : sampledSurface section3
# Faces         : 10
# Area          : 6.45160004e-02
# Scale factor  : 1.00000000e+00
# Time          	areaAverage(alpha.liquid)	areaAverage(p)
0.0031113797    	0.00000000e+00	8.55607419e+06
0.36329334      	2.87188631e-86	8.55584009e+06
0.72066189      	3.80204855e-76	8.55555631e+06
1.0798016       	5.03902112e-66	8.55567798e+06
1.3989534       	3.29158021e-57	8.55560592e+06
1.7042669       	6.27455731e-49	8.55563271e+06
2.0054099       	6.15949443e-41	8.55565086e+06
2.3105285       	4.44779917e-33	8.55565581e+06
2.6222961       	2.13754539e-25	8.55566397e+06
2.9409577       	2.16122051e-17	8.55565446e+06
3.3149079       	2.91146091e-11	8.55565939e+06
3.8034098       	2.37459408e-06	8.55566043e+06
4.3920279       	2.67380558e-03	8.55563305e+06
5.0183699       	8.89978804e-04	8.55561610e+06
5.5596276       	3.13913166e-01	8.55627582e+06
5.9155984       	6.89359774e-01	8.55666838e+06
6.2102353       	6.00585063e-01	8.55661716e+06
6.4926321       	6.47060602e-01	8.55656855e+06
6.7476079       	5.42749130e-01	8.55655304e+06
7.0082041       	6.14313108e-01	8.55476311e+06
7.2625014       	6.59437609e-01	8.55605912e+06
7.5019608       	7.61391408e-01	8.55548235e+06
7.7362922       	7.40877761e-01	8.55732668e+06
7.9856821       	7.01136428e-01	8.55766185e+06
8.2618065       	6.60820312e-01	8.55692856e+06
8.5415587       	7.95385407e-01	8.55687031e+06
8.7844041       	7.03666996e-01	8.55769443e+06
9.0300112       	7.54105176e-01	8.55865605e+06
9.2673174       	7.90894466e-01	8.55911328e+06
9.5113947       	7.78311317e-01	8.55879632e+06
9.7491943       	6.49868058e-01	8.55679056e+06
9.9847677       	6.53605331e-01	8.55772145e+06
10.228487       	7.32363773e-01	8.55770500e+06
10.473112       	6.64702749e-01	8.55784713e+06
10.715338       	7.61676513e-01	8.55787711e+06
10.962178       	8.03634769e-01	8.55860978e+06
11.203224       	9.27758606e-01	8.55961441e+06
11.45013        	6.06482052e-01	8.55754148e+06
11.69024        	7.67980559e-01	8.55869294e+06
11.929091       	7.34935790e-01	8.55823278e+06
12.170959       	7.25065440e-01	8.55736550e+06
12.418552       	7.48561250e-01	8.55850735e+06
12.659236       	8.70280901e-01	8.55761853e+06
12.898039       	7.45978917e-01	8.55800359e+06
13.141231       	7.54340911e-01	8.55839766e+06
13.383018       	8.63481039e-01	8.55795303e+06
13.624803       	6.75615556e-01	8.55660355e+06
13.870181       	8.68939784e-01	8.55861044e+06
14.117744       	7.80630693e-01	8.55850964e+06
14.35358        	7.63426683e-01	8.55951763e+06
14.59921        	6.96919511e-01	8.55939039e+06
14.837704       	7.30315905e-01	8.55652686e+06
15.080992       	7.47805928e-01	8.55761442e+06
15.326917       	8.24841496e-01	8.55907514e+06
15.566603       	8.54259739e-01	8.55882790e+06
15.803277       	6.69054802e-01	8.55641311e+06
16.048048       	7.70604958e-01	8.55793825e+06
16.292131       	7.17596985e-01	8.55854694e+06
16.532814       	6.15335911e-01	8.55805965e+06
16.776131       	8.42560067e-01	8.55832995e+06
17.020316       	5.31168668e-01	8.55642714e+06
17.255744       	7.10518323e-01	8.55792452e+06
17.496635       	6.87825457e-01	8.55823895e+06
17.734876       	7.35850085e-01	8.55774094e+06
17.979185       	7.93542924e-01	8.55777393e+06
18.222654       	6.08723826e-01	8.55729375e+06
18.470661       	8.20391908e-01	8.55847190e+06
18.716408       	6.78199973e-01	8.55778113e+06
18.965916       	8.76921906e-01	8.55755523e+06
19.21266        	7.67503208e-01	8.55783498e+06
19.453065       	7.32051138e-01	8.55768105e+06
19.692758       	7.09116217e-01	8.55749058e+06
19.941695       	7.38698300e-01	8.55688143e+06
20.188445       	7.00019244e-01	8.55829089e+06
20.433951       	7.62591100e-01	8.55818849e+06
20.679175       	8.72310389e-01	8.55986069e+06
20.929055       	9.31937615e-01	8.55923317e+06
21.176505       	7.53162813e-01	8.55831035e+06
21.415867       	8.57233050e-01	8.55825608e+06
21.663796       	6.86260337e-01	8.55852906e+06
21.909817       	7.34453680e-01	8.55818983e+06
22.156909       	8.38441225e-01	8.55956707e+06
22.403709       	7.78831832e-01	8.55740413e+06

`,
  // DPMFoam with the particleErosion cloud function: the generated erosion case run in OpenFOAM v2412 (first rows of each file)
  erosion: `### file: postProcessing/cloudInfo/0/kinematicCloud.dat
# Cloud information
# Time          	nParcels        	mass            	Dmax            	D10             	D32             
0.021681352     	434	1.25470765e-05	1.50000000e-04	1.50000000e-04	1.50000000e-04
0.071839137     	1435	4.15735673e-05	1.50000000e-04	1.50000000e-04	1.50000000e-04
0.13888889      	2776	8.03755000e-05	1.50000000e-04	1.50000000e-04	1.50000000e-04
0.20833333      	4165	1.20563250e-04	1.50000000e-04	1.50000000e-04	1.50000000e-04
0.27777778      	5554	1.60751000e-04	1.50000000e-04	1.50000000e-04	1.50000000e-04
0.34722222      	6943	2.00938750e-04	1.50000000e-04	1.50000000e-04	1.50000000e-04
0.41666667      	8332	2.41126500e-04	1.50000000e-04	1.50000000e-04	1.50000000e-04
0.48611111      	9720	2.81314250e-04	1.50000000e-04	1.50000000e-04	1.50000000e-04
0.55555556      	11108	3.21502000e-04	1.50000000e-04	1.50000000e-04	1.50000000e-04
0.625           	12498	3.61689750e-04	1.50000000e-04	1.50000000e-04	1.50000000e-04
0.69444444      	13886	4.01877500e-04	1.50000000e-04	1.50000000e-04	1.50000000e-04
0.76388889      	15251	4.41341427e-04	1.50000000e-04	1.50000000e-04	1.50000000e-04
0.83333333      	16334	4.72695436e-04	1.50000000e-04	1.50000000e-04	1.50000000e-04
0.90277778      	16998	4.91891990e-04	1.50000000e-04	1.50000000e-04	1.50000000e-04
0.97222222      	17470	5.05509102e-04	1.50000000e-04	1.50000000e-04	1.50000000e-04
1.0416667       	17809	5.15280637e-04	1.50000000e-04	1.50000000e-04	1.50000000e-04
1.1111111       	18101	5.23728332e-04	1.50000000e-04	1.50000000e-04	1.50000000e-04
1.1805556       	18361	5.31280095e-04	1.50000000e-04	1.50000000e-04	1.50000000e-04
1.25            	18602	5.38247506e-04	1.50000000e-04	1.50000000e-04	1.50000000e-04
1.3194444       	18936	5.47927951e-04	1.50000000e-04	1.50000000e-04	1.50000000e-04
1.3888889       	19150	5.54083366e-04	1.50000000e-04	1.50000000e-04	1.50000000e-04
1.4583333       	19395	5.61160445e-04	1.50000000e-04	1.50000000e-04	1.50000000e-04
1.5277778       	19628	5.67903558e-04	1.50000000e-04	1.50000000e-04	1.50000000e-04
1.5972222       	19883	5.75281268e-04	1.50000000e-04	1.50000000e-04	1.50000000e-04
1.6666667       	20085	5.81151267e-04	1.50000000e-04	1.50000000e-04	1.50000000e-04
1.7361111       	20262	5.86305790e-04	1.50000000e-04	1.50000000e-04	1.50000000e-04
### file: postProcessing/erosionMax/0/surfaceFieldValue.dat
# Region type : patch wall
# Faces         : 408
# Area          : 1.27088242e+00
# Scale factor  : 1.00000000e+00
# Time          	max(kinematicCloudQ)
0.021681352     	0.00000000e+00
0.071839137     	2.54987988e-19
0.13888889      	9.85292466e-19
0.20833333      	1.69598875e-18
0.27777778      	2.80971896e-18
0.34722222      	3.65080596e-18
0.41666667      	4.79243465e-18
0.48611111      	5.13666934e-18
0.55555556      	5.68727780e-18
0.625           	6.65713754e-18
0.69444444      	7.55449966e-18
0.76388889      	8.96098025e-18
0.83333333      	1.00684744e-17
0.90277778      	1.07012870e-17
0.97222222      	1.09776476e-17
1.0416667       	1.13078918e-17
1.1111111       	1.16744864e-17
1.1805556       	1.18201433e-17
1.25            	1.21014737e-17
1.3194444       	1.24157514e-17
1.3888889       	1.31465579e-17
1.4583333       	1.40649640e-17
1.5277778       	1.42813568e-17
### file: postProcessing/erosionTotal/0/surfaceFieldValue.dat
# Region type : patch wall
# Faces         : 408
# Area          : 1.27088242e+00
# Scale factor  : 1.00000000e+00
# Time          	sum(kinematicCloudQ)
0.021681352     	0.00000000e+00
0.071839137     	4.18896715e-19
0.13888889      	1.98269123e-18
0.20833333      	5.23149221e-18
0.27777778      	1.12645991e-17
0.34722222      	1.62673033e-17
0.41666667      	2.25462204e-17
0.48611111      	3.37404169e-17
0.55555556      	4.21521974e-17
0.625           	5.48670814e-17
0.69444444      	7.10690342e-17
0.76388889      	8.60221617e-17
0.83333333      	1.04734373e-16
0.90277778      	1.24110848e-16
0.97222222      	1.51363574e-16
1.0416667       	1.83286097e-16
1.1111111       	2.19165557e-16
1.1805556       	2.57782347e-16
1.25            	3.01605029e-16
1.3194444       	3.45662460e-16
1.3888889       	3.95123396e-16
1.4583333       	4.49831403e-16
1.5277778       	5.07888753e-16
### file: postProcessing/inletPressure/0/surfaceFieldValue.dat
# Region type : patch inlet
# Faces         : 45
# Area          : 4.83870000e-02
# Scale factor  : 1.00000000e+00
# Time          	areaAverage(pStatic)
0.021681352     	4.75354568e+03
0.071839137     	4.46171065e+03
0.13888889      	4.43294868e+03
0.20833333      	4.43417073e+03
0.27777778      	4.43720924e+03
0.34722222      	4.44008043e+03
0.41666667      	4.44321920e+03
0.48611111      	4.44681994e+03
0.55555556      	4.45080279e+03
0.625           	4.45478582e+03
0.69444444      	4.45820704e+03
0.76388889      	4.46093487e+03
0.83333333      	4.46298729e+03
0.90277778      	4.46446231e+03
0.97222222      	4.46545618e+03
1.0416667       	4.46609846e+03
1.1111111       	4.46648667e+03
1.1805556       	4.46671105e+03
1.25            	4.46687688e+03
1.3194444       	4.46696339e+03
1.3888889       	4.46702032e+03
1.4583333       	4.46705190e+03
1.5277778       	4.46708135e+03
### file: postProcessing/outletPressure/0/surfaceFieldValue.dat
# Region type : patch outlet
# Faces         : 45
# Area          : 4.83870000e-02
# Scale factor  : 1.00000000e+00
# Time          	areaAverage(pStatic)
0.021681352     	0.00000000e+00
0.071839137     	0.00000000e+00
0.13888889      	0.00000000e+00
0.20833333      	0.00000000e+00
0.27777778      	0.00000000e+00
0.34722222      	0.00000000e+00
0.41666667      	0.00000000e+00
0.48611111      	0.00000000e+00
0.55555556      	0.00000000e+00
0.625           	0.00000000e+00
0.69444444      	0.00000000e+00
0.76388889      	0.00000000e+00
0.83333333      	0.00000000e+00
0.90277778      	0.00000000e+00
0.97222222      	0.00000000e+00
1.0416667       	0.00000000e+00
1.1111111       	0.00000000e+00
1.1805556       	0.00000000e+00
1.25            	0.00000000e+00
1.3194444       	0.00000000e+00
1.3888889       	0.00000000e+00
1.4583333       	0.00000000e+00
1.5277778       	0.00000000e+00
`,
  // CalculiX 2.23, pipe_modal.dat of the generated shell span (6 m, fixed ends)
  modalDat: `
                        S T E P       1


     E I G E N V A L U E   O U T P U T

 MODE NO    EIGENVALUE                       FREQUENCY   
                                     REAL PART            IMAGINARY PART
                           (RAD/TIME)      (CYCLES/TIME     (RAD/TIME)

      1   0.4727589E+05   0.2174302E+03   0.3460509E+02   0.0000000E+00
      2   0.4727589E+05   0.2174302E+03   0.3460509E+02   0.0000000E+00
      3   0.3339410E+06   0.5778763E+03   0.9197187E+02   0.0000000E+00
      4   0.3339410E+06   0.5778763E+03   0.9197187E+02   0.0000000E+00
      5   0.1175470E+07   0.1084191E+04   0.1725544E+03   0.0000000E+00
      6   0.1175470E+07   0.1084191E+04   0.1725544E+03   0.0000000E+00
      7   0.1498591E+07   0.1224170E+04   0.1948327E+03   0.0000000E+00
      8   0.2908307E+07   0.1705376E+04   0.2714190E+03   0.0000000E+00
      9   0.2908307E+07   0.1705376E+04   0.2714190E+03   0.0000000E+00
     10   0.3902844E+07   0.1975562E+04   0.3144204E+03   0.0000000E+00

     P A R T I C I P A T I O N   F A C T O R S
`,
  // CalculiX 2.23, pipe_static.dat of the generated solid bend (first displacement block and reaction totals)
  staticDat: `
                        S T E P       1


                                INCREMENT     1


 displacements (vx,vy,vz) for set BEND and time  0.1000000E+01

      1093 -1.526780E-04 -3.685145E-05  1.526780E-04
      1094 -1.527978E-04 -3.636286E-05  1.527978E-04
      1095 -1.529185E-04 -3.593192E-05  1.529185E-04
      1096 -1.485093E-04 -3.563541E-05  1.485093E-04
      1097 -1.481872E-04 -3.447938E-05  1.481872E-04

 total force (fx,fy,fz) for set END0 and time  0.1000000E+01
        4.124111E+05  3.782770E-08 -1.236894E+03
 total force (fx,fy,fz) for set END0 and time  0.2000000E+01
        4.194147E+05  1.919252E-07  6.640295E+03
`,
  // CalculiX 2.23, pipe_static.frd of the generated shell span (header, and the first four nodes of every block)
  frd: `    1C
    1UUSER                                                              
    1UDATE              08.october.2026                                 
    1UTIME              21:11:35                                        
    1UHOST                                                              
    1UPGM               CalculiX                                        
    1UVERSION           Version 2.23                             
    1UCOMPILETIME       Mon Sep 28 09:16:29 2026                    
    1UDIR                                                               
    1UDBN                                                               
    1UMAT    1STEEL                                                     
    2C                          2076                                     1
 -1       889 0.00000E+00 1.27000E-01 0.00000E+00
 -1       890 0.00000E+00 1.34950E-01 0.00000E+00
 -1       891 0.00000E+00 1.42900E-01 0.00000E+00
 -1       892 1.70511E-18 1.22673E-01 3.28700E-02
 -3
    3C                           288                                     1
 -1         1    4    0    1
 -1         2    4    0    1
 -1         3    4    0    1
 -1         4    4    0    1
 -3
    1PSTEP                         1           1           1          
  100CL  101 1.000000000        2076                     0    1           1
 -4  DISP        4    1
 -5  D1          1    2    1    0
 -5  D2          1    2    2    0
 -5  D3          1    2    3    0
 -5  ALL         1    2    0    0    1ALL
 -1       889 0.00000E+00-5.23154E-07 0.00000E+00
 -1       890 0.00000E+00 0.00000E+00 0.00000E+00
 -1       891 0.00000E+00 5.23154E-07 0.00000E+00
 -1       892 0.00000E+00-1.23189E-07-3.30084E-08
 -3
    1PSTEP                         2           1           1          
  100CL  101 1.000000000        2076                     0    1           1
 -4  STRESS      6    1
 -5  SXX         1    4    1    1
 -5  SYY         1    4    2    2
 -5  SZZ         1    4    3    3
 -5  SXY         1    4    1    2
 -5  SYZ         1    4    2    3
 -5  SZX         1    4    3    1
 -1       889-2.81312E+07-4.30262E+06-2.87820E+07-5.31798E+06-1.17859E-06 9.00832E-06
 -1       890-1.35035E+07 3.11917E+06-2.11440E+07-4.70895E+06-7.77305E-07 8.67335E-06
 -1       891 1.12424E+06 1.05410E+07-1.35061E+07-4.09992E+06-3.76021E-07 8.33836E-06
 -1       892-2.81312E+07-6.00871E+06-2.70759E+07-4.95632E+06 6.08157E+06-1.32804E+06
 -3
    1PSTEP                         3           1           1          
  100CL  101 1.000000000        2076                     0    1           1
 -4  ERROR       1    1
 -5  STR(%)      1    1    0    0
 -1       889 2.95178E+01
 -1       890 2.95178E+01
 -1       891 2.95178E+01
 -1       892 2.95178E+01
 -3
    1PSTEP                         4           1           2          
  100CL  102 2.000000000        2076                     0    2           1
 -4  DISP        4    1
 -5  D1          1    2    1    0
 -5  D2          1    2    2    0
 -5  D3          1    2    3    0
 -5  ALL         1    2    0    0    1ALL
 -1       889 0.00000E+00-8.21732E-06 0.00000E+00
 -1       890 0.00000E+00 0.00000E+00 0.00000E+00
 -1       891 0.00000E+00 8.21732E-06 0.00000E+00
 -1       892 0.00000E+00-9.00799E-06-2.41368E-06
 -3
    1PSTEP                         5           1           2          
  100CL  102 2.000000000        2076                     0    2           1
 -4  STRESS      6    1
 -5  SXX         1    4    1    1
 -5  SYY         1    4    2    2
 -5  SZZ         1    4    3    3
 -5  SXY         1    4    1    2
 -5  SYZ         1    4    2    3
 -5  SZX         1    4    3    1
 -1       889-1.11373E+08 2.27770E+07-1.29921E+08 1.52296E+07 5.86861E+04-2.69510E+05
 -1       890-1.45141E+08 9.44973E+06-1.31518E+08 1.34772E+07 3.98795E+04-2.17661E+05
 -1       891-1.78909E+08-3.87750E+06-1.33116E+08 1.17247E+07 2.10729E+04-1.65812E+05
 -1       892-1.10954E+08 1.24701E+06-1.08369E+08 1.42392E+07 3.16410E+07 3.55580E+06
 -3
    1PSTEP                         6           1           2          
  100CL  102 2.000000000        2076                     0    2           1
 -4  ERROR       1    1
 -5  STR(%)      1    1    0    0
 -1       889 1.00874E+01
 -1       890 1.00874E+01
 -1       891 1.00874E+01
 -1       892 1.00935E+01
 -3
 9999

`,
  // first rows of gerg2008_table.csv written by the generated CoolProp script
  gergHead: `P,T,wG,rhoG,zG,cpG,jtG,wSound,rhoO,cpO,source
1,-30,,1.0644305,0.99414975,1811.0465,9.4742481e-06,346.01004,,,GERG-2008 (CoolProp 6.7.0)
1,-17.5,,1.0463175,0.99461832,1818.7298,8.924863e-06,347.19723,,,GERG-2008 (CoolProp 6.7.0)
1,-5,,1.0301813,0.99503247,1833.103,8.393846e-06,348.07083,,,GERG-2008 (CoolProp 6.7.0)
`,
};

// ---------------------------------------------------------------------------------------------- hand-off table
section('hand-off table');
// items the suite engineers asked to be handed off (the request files at the time of writing)
const REQUESTED = { flow: ['LES', 'DES', 'IDDES', 'DNS where computationally feasible', '1-D transient flow + 3-D CFD', 'LES + VOF'], integ: ['Navier–Stokes + structural dynamics', 'CFD + FEA'], solids: ['Population balance + CFD'] };
// items named in the brief as beyond a browser
const MUST = { flow: ['LES', 'DES', 'IDDES', 'DNS where computationally feasible', 'LES + VOF', '1-D transient flow + 3-D CFD'], integ: ['CFD + FEA', 'Navier–Stokes + structural dynamics', 'Erosion + particle CFD'], solids: ['Population balance + CFD'], pvt: ['GERG-type multiparameter EOS'] };
ok(SUITE_IDS.every((id) => Array.isArray(B.HANDOFF[id])), 'HANDOFF has an array for every suite');
for (const id of SUITE_IDS) {
  const items = itemsOf(id), norms = items.map(B.normItem), listed = B.HANDOFF_ITEMS(id);
  ok(new Set(B.HANDOFF[id].map((e) => e.match)).size === B.HANDOFF[id].length, `${id}: no duplicate entries`);
  for (const e of B.HANDOFF[id]) {
    ok(items.includes(e.item), `${id}: "${e.item}" is an exact catalogue item of this suite`);
    ok(e.match === B.normItem(e.item) && norms.filter((n) => n === e.match).length === 1, `${id}: "${e.item}" matches exactly one catalogue item`);
    ok(typeof B[e.generator] === 'function', `${id}: generator ${e.generator} of "${e.item}" exists`);
    ok(typeof e.solver === 'string' && e.solver.length > 5 && typeof e.what === 'string' && e.what.length > 30, `${id}: "${e.item}" names a solver and says what it resolves`);
    ok(B.handoffFor(id, e.item) === e, `${id}: handoffFor finds "${e.item}"`);
  }
  // nothing else in the catalogue of the suite matches (no short-fragment accidents such as "les" or "des")
  for (const it of items) ok((B.handoffFor(id, it) !== null) === listed.includes(it), `${id}: "${it}" ${listed.includes(it) ? 'is' : 'is not'} handed off`);
  for (const it of [...(REQUESTED[id] || []), ...(MUST[id] || [])]) ok(B.handoffFor(id, it) !== null, `${id}: requested item "${it}" has an entry`);
}
ok(B.handoffFor('flow', 'Slug-capturing model') === null && B.handoffFor('flow', 'Mandhane-type flow maps') === null && B.handoffFor('flow', 'Slug-length correlations') === null, 'items containing "les"/"des" inside other words are not matched');
ok(B.handoffFor('nosuchsuite', 'LES') === null && B.handoffFor('flow', '') === null, 'unknown suite or empty name gives null');
ok(B.handoffFor('flow', 'les') !== null && B.handoffFor('integ', 'Navier-Stokes + structural dynamics') !== null, 'matching ignores case and dash style');
{ // request files, when the folder is there
  const dir = process.env.HANDOFF_REQUESTS || '';
  let n = 0;
  if (dir && fs.existsSync(dir)) for (const f of fs.readdirSync(dir).filter((x) => x.endsWith('.json'))) {
    const id = f.replace(/\.json$/, '');
    let list = []; try { list = JSON.parse(fs.readFileSync(path.join(dir, f), 'utf8')); } catch { ok(false, `request file ${f} is not valid JSON`); }
    for (const r of Array.isArray(list) ? list : []) { n++; ok(itemsOf(id).includes(r.item), `request ${f}: "${r.item}" is a catalogue item`); ok(B.handoffFor(id, r.item) !== null, `request ${f}: "${r.item}" has a hand-off entry`); }
  }
  console.log(`   request files: ${n} requested item${n === 1 ? '' : 's'} checked`);
}

// ---------------------------------------------------------------------------------------------- zip
section('zip writer');
const enc = new TextEncoder();
ok(B.crc32(enc.encode('123456789')) === 0xcbf43926, 'CRC-32 check value');
ok(B.crc32(new Uint8Array(0)) === 0, 'CRC-32 of nothing');
{
  const files = [{ path: 'Allrun', text: '#!/bin/sh\necho hi\n' }, { path: 'system/controlDict', text: 'application interFoam;\n' }, { path: 'data/naïve–name.txt', text: 'ü — ✓\n' }, { path: 'empty', text: '' }, { path: 'bin.dat', bytes: new Uint8Array([0, 255, 1, 254]) }];
  const z = B.zipStore(files, { root: 'case' }), dv = new DataView(z.buffer);
  ok(dv.getUint32(0, true) === 0x04034b50, 'starts with a local file header');
  ok(dv.getUint32(z.length - 22, true) === 0x06054b50 && dv.getUint16(z.length - 12, true) === files.length, 'ends with the end record listing every entry');
  const cdOff = dv.getUint32(z.length - 6, true), cdSize = dv.getUint32(z.length - 10, true);
  ok(cdOff + cdSize + 22 === z.length && dv.getUint32(cdOff, true) === 0x02014b50, 'central directory offset and size are consistent');
  const back = B.unzipStored(z);
  ok(back.length === files.length && back.every((e) => e.crcOk && e.method === 0), 'every entry is stored and its CRC verifies');
  ok(back.map((e) => e.path).join('|') === files.map((f) => 'case/' + f.path).join('|'), 'paths survive (UTF-8, root folder)');
  ok(back[2].text === 'ü — ✓\n' && back[0].text === files[0].text && back[4].bytes.join() === '0,255,1,254', 'contents survive');
  ok((back[0].mode & 0o111) !== 0 && (back[1].mode & 0o111) === 0, 'run scripts carry the executable bit, dictionaries do not');
  ok(B.zipStore(files).join() === B.zipStore(files).join(), 'archive is reproducible');
  let threw = 0; for (const bad of [[{ path: '../x', text: '' }], [{ path: 'a', text: '' }, { path: 'a', text: '' }]]) { try { B.zipStore(bad); } catch { threw++; } }
  ok(threw === 2, 'unsafe and duplicate paths are refused');
  const cut = z.slice(0, z.length - 30); let t2 = false; try { B.unzipStored(cut); } catch { t2 = true; } ok(t2, 'a truncated archive is refused');
  // an independent reader, when the system has one
  try { const dir = fs.mkdtempSync(path.join(os.tmpdir(), 'hsbridge-')), f = path.join(dir, 'c.zip'); fs.writeFileSync(f, z); const out = execFileSync('unzip', ['-t', f], { encoding: 'utf8', stdio: ['ignore', 'pipe', 'ignore'] }); ok(/No errors detected/.test(out), 'unzip -t accepts the archive'); fs.rmSync(dir, { recursive: true, force: true }); } catch (e) { if (e.code !== 'ENOENT') ok(false, 'unzip -t rejected the archive: ' + String(e.message).slice(0, 200)); }
}

// ---------------------------------------------------------------------------------------------- meshes
section('block meshes');
const D = 0.254;
{
  const cases = { straight: [{ type: 'straight', L: 5 }], bend90: [{ type: 'straight', L: 2 }, { type: 'bend', R: 0.4, angleDeg: 90, toward: 'up' }, { type: 'straight', L: 2 }], bend135left: [{ type: 'straight', L: 1 }, { type: 'bend', R: 1.2, angleDeg: 135, toward: 'left' }, { type: 'straight', L: 1 }], jumper: [{ type: 'straight', L: 3 }, { type: 'bend', R: 0.8, angleDeg: 90, toward: 'down' }, { type: 'straight', L: 6 }, { type: 'bend', R: 0.8, angleDeg: 90, toward: 'down' }, { type: 'straight', L: 3 }] };
  for (const [name, legs] of Object.entries(cases)) {
    const L = legs.reduce((s, l) => s + (l.type === 'bend' ? (l.R * l.angleDeg * Math.PI) / 180 : l.L), 0);
    for (const m of [B.ogridMesh({ D, legs, nc: 6, nr: 5, wallRatio: 0.3, dz: 0.05 }), B.slotMesh({ D, legs, nD: 12, dz: 0.05 })]) {
      const c = B.checkBlockMesh(m), exact = m.type === 'ogrid' ? (Math.PI / 4) * D * D * L : D * m.width * L;
      ok(c.ok, `${name} ${m.type}: topology (${c.issues.slice(0, 2).join('; ')})`);
      ok(c.minJacobian > 0, `${name} ${m.type}: positive Jacobian at every block corner`);
      near(m.length, L, 1e-12, `${name} ${m.type}: developed length`);
      near(B.blockMeshVolume(m), exact, 2e-4, `${name} ${m.type}: volume equals section area × length`);
      // closed O-grid: 5 blocks per segment, every station ring closed, Euler count of the block complex
      if (m.type === 'ogrid') { const seg = m.stations.length - 1; ok(m.blocks.length === 5 * seg && m.vertices.length === 8 * (seg + 1) && c.internalFaces === 8 * seg + 5 * (seg - 1) && c.boundaryFaces === 4 * seg + 10, `${name}: closed O-grid block structure`); ok(m.nCells === m.blocks.reduce((s, b) => s + b.n[0] * b.n[1] * b.n[2], 0), `${name}: cell count`); }
      const txt = B.blockMeshDict(m), nHex = (txt.match(/^\s+hex \(/gm) || []).length, nArc = (txt.match(/^\s+arc /gm) || []).length;
      ok(nHex === m.blocks.length && nArc === m.edges.length && /FoamFile/.test(txt) && /mergePatchPairs/.test(txt), `${name} ${m.type}: blockMeshDict lists every block and edge`);
      ok((txt.match(/[({]/g) || []).length === (txt.match(/[)}]/g) || []).length, `${name} ${m.type}: brackets balance in blockMeshDict`);
    }
  }
  const end = B.sweepStations(cases.jumper).pop();
  ok(Math.abs(end.c[0]) < 1e-9 && Math.abs(end.c[1]) < 1e-9 && Math.abs(end.c[2] + 6 + 1.6) < 1e-9 && Math.abs(end.t[0] + 1) < 1e-12, 'jumper path ends level with its start, pointing back down');
  near(B.levelForHoldup(0.5), 0.5, 1e-9, 'half-full pipe level'); ok(B.levelForHoldup(0.1) < 0.2 && B.levelForHoldup(0.9) > 0.8, 'stratified level is monotone');
}
section('finite-element meshes');
for (const element of ['shell', 'solid', 'solid8']) for (const legs of [[{ type: 'straight', L: 4 }], [{ type: 'straight', L: 1.5 }, { type: 'bend', R: 1.27, angleDeg: 90, toward: 'up' }, { type: 'straight', L: 1.5 }]]) {
  const m = B.pipeFeMesh({ D, wt: 0.0159, legs, element, nC: 16, elemLen: 0.15 }), c = B.checkFeMesh(m), tag = `${element} ${legs.length > 1 ? 'bend' : 'span'}`;
  ok(c.ok, `${tag}: every element references existing nodes, positive Jacobians (${c.issues.slice(0, 2).join('; ')})`);
  near(c.volume, c.exactVolume, element === 'solid8' ? 0.03 : 5e-4, `${tag}: wall volume π (ro² − ri²) L`);
  ok(m.elements.every((e) => e.length === (element === 'shell' ? 8 : element === 'solid' ? 20 : 8)), `${tag}: nodes per element`);
  ok(m.sets.END0.length === m.sets.END1.length && m.sets.END0.length > 0 && m.sets.MID.length > 0 && m.sets.WETTED.length > 0 && (legs.length > 1) === (m.sets.BEND.length > 0), `${tag}: node sets`);
  const inp = B.calculixMeshText(m), nNode = (inp.slice(inp.indexOf('*NODE'), inp.indexOf('*ELEMENT')).match(/^\d+, -?[\d.e-]+, -?[\d.e-]+, -?[\d.e-]+$/gm) || []).length;
  ok(nNode === m.nodes.length && inp.includes(`TYPE=${m.type}`), `${tag}: mesh include file lists every node`);
}
near(B.beamFrequency(6, 207e9 * 1.1e-4, 150, 'pinned'), (Math.PI / 2 / 36) * Math.sqrt((207e9 * 1.1e-4) / 150), 1e-12, 'beam frequency, pinned–pinned');

// ---------------------------------------------------------------------------------------------- generators
section('generators on the default case');
const snap = B.bridgeCase({});
ok(snap.source === 'kernel estimate' && snap.D === 0.254 && snap.holdup > 0 && snap.holdup < 1 && snap.rhoL > snap.rhoG && snap.slug.freq > 0 && snap.slug.forceN > 0, 'case snapshot from the kernel estimate (nothing solved yet)');
ok(Math.abs(B.bridgeCase({}, { x: 19300 }).angleDeg) > 30 && Math.abs(B.bridgeCase({}, { x: 5000 }).angleDeg) < 5, 'local inclination follows the profile');
const grid = B.propertyGrid(DEFAULT_FLUID, { P: [5, 50, 200], T: [10, 60] });
const EXPECT = {
  openfoamPipeCase: [{}, ['0.orig/U', '0.orig/p_rgh', '0.orig/alpha.liquid', '0.orig/k', '0.orig/omega', '0.orig/nut', 'constant/transportProperties', 'constant/turbulenceProperties', 'constant/g', 'system/blockMeshDict', 'system/controlDict', 'system/fvSchemes', 'system/fvSolution', 'system/decomposeParDict', 'system/setFieldsDict', 'system/topoSetDict', 'system/createPatchDict', 'Allrun', 'Allclean', 'collect.sh', 'README.md']],
  openfoamCoupledCase: [{}, ['0.orig/U', 'constant/inletSeries.csv', 'system/controlDict', 'Allrun', 'README.md']],
  openfoamRiserSlugCase: [{}, ['0.orig/U', '0.orig/p', '0.orig/T', '0.orig/alphat', 'constant/thermophysicalProperties', 'constant/thermophysicalProperties.gas', 'constant/thermophysicalProperties.liquid', 'system/blockMeshDict', 'Allrun', 'README.md']],
  openfoamErosionCase: [{}, ['0.orig/U.carrier', '0.orig/p', 'constant/kinematicCloudProperties', 'constant/transportProperties', 'system/blockMeshDict', 'system/controlDict', 'Allrun', 'collect.sh', 'README.md']],
  openfoamEulerCase: [{}, ['constant/phaseProperties', 'system/blockMeshDict', 'system/controlDict', 'system/fvSchemes', 'system/fvSolution', 'Allrun', 'collect.sh', 'README.md']],
  openfoamPopulationCase: [{}, ['constant/phaseProperties', 'system/blockMeshDict', 'system/controlDict', 'Allrun', 'collect.sh', 'README.md']],
  calculixPipeCase: [{}, ['mesh.inc', 'pipe_static.inp', 'pipe_modal.inp', 'pipe_dynamic.inp', 'run.sh', 'README.md']],
  fsiCase: [{}, ['precice-config.xml', 'fluid-openfoam/system/preciceDict', 'fluid-openfoam/system/blockMeshDict', 'fluid-openfoam/system/controlDict', 'fluid-openfoam/constant/dynamicMeshDict', 'fluid-openfoam/0/pointDisplacement', 'fluid-openfoam/0/U', 'fluid-openfoam/0/p', 'solid-calculix/pipe.inp', 'solid-calculix/mesh.inc', 'solid-calculix/config.yml', 'Allrun', 'README.md']],
  coolpropScript: [{ grid }, ['gerg2008_coolprop.py', 'README.md']],
  teqpScript: [{ grid }, ['pcsaft_teqp.py', 'README.md']],
};
const FORBIDDEN = /\b(NaN|undefined|Infinity)\b|\[object Object\]/, built = {};
for (const [gen, [extra, expected]] of Object.entries(EXPECT)) {
  const make = () => B[gen]({ case: snap, ...extra }), r = make(); built[gen] = r;
  const paths = r.files.map((f) => f.path);
  ok(expected.every((p) => paths.includes(p)), `${gen}: expected files (missing: ${expected.filter((p) => !paths.includes(p)).join(', ') || 'none'})`);
  ok(new Set(paths).size === paths.length && r.files.every((f) => typeof f.text === 'string'), `${gen}: unique paths, text contents`);
  const bad = r.files.filter((f) => FORBIDDEN.test(f.text)).map((f) => f.path);
  ok(!bad.length, `${gen}: no NaN / undefined / Infinity in any file (${bad.join(', ')})`);
  ok(typeof r.readme === 'string' && r.readme.length > 400 && typeof r.summary === 'string' && r.summary.length > 20 && Array.isArray(r.commands) && r.commands.length >= 2 && typeof r.name === 'string', `${gen}: read-me, summary, commands, name`);
  ok(r.files.find((f) => f.path === 'README.md')?.text === r.readme, `${gen}: README.md is the returned read-me`);
  ok(JSON.stringify(make().files) === JSON.stringify(r.files), `${gen}: deterministic`);
  ok(r.files.every((f) => f.path.includes('README') || !/^#!/.test(f.text) || /^#!\/(bin\/sh|usr\/bin\/env python3)\n/.test(f.text)), `${gen}: scripts start with a portable interpreter line`);
  for (const f of r.files) if (/(Dict|Properties|^0\.orig\/|\/0\/|fvS)/.test(f.path) && /FoamFile/.test(f.text)) { const t = f.text.replace(/\/\/.*$/gm, '').replace(/"[^"\n]*"/g, '""'); if ((t.match(/\{/g) || []).length !== (t.match(/\}/g) || []).length || (t.match(/\(/g) || []).length !== (t.match(/\)/g) || []).length) ok(false, `${gen}: unbalanced brackets in ${f.path}`); }
  const z = B.unzipStored(B.zipStore(r.files, { root: r.name })); ok(z.length === r.files.length && z.every((e) => e.crcOk), `${gen}: zips and reads back`);
  ok(typeof B[gen]().summary === 'string', `${gen}: works without any input (reference numbers)`);
}
ok(Object.values(B.HANDOFF).flat().every((e) => EXPECT[e.generator]), 'every generator named in the hand-off table is exercised here');
for (const e of Object.values(B.HANDOFF).flat()) { let r = null; try { r = B[e.generator]({ case: snap, grid, ...e.options }); } catch (err) { ok(false, `"${e.item}": ${err.message}`); } ok(r && r.files.length > 1, `"${e.item}": its options give a case`); }
{ // turbulence models, interface methods, geometries
  for (const key of Object.keys(B.TURBULENCE)) { const r = B.openfoamPipeCase({ case: snap, turbulence: key, cellsPerDiameter: 12, lengthD: 5 }), tp = r.files.find((f) => f.path === 'constant/turbulenceProperties').text, tm = B.TURBULENCE[key]; ok(tm.type === 'laminar' ? /simulationType\s+laminar/.test(tp) : tp.includes(tm.type + 'Model') && tp.includes(key) && tm.fields.every((f) => r.files.some((x) => x.path === '0.orig/' + f)), `turbulence ${key}: model and its fields are written`); ok(!/IDDES$/.test(key) || /IDDESDelta/.test(tp), `turbulence ${key}: IDDES delta`); }
  const iso = B.openfoamPipeCase({ case: snap, interface: 'isoAdvector' }), rdf = B.openfoamPipeCase({ case: snap, interface: 'plicRDF' });
  ok(/application\s+interIsoFoam/.test(iso.files.find((f) => f.path === 'system/controlDict').text) && /isoAlpha/.test(iso.files.find((f) => f.path === 'system/fvSolution').text) && /plicRDF/.test(rdf.files.find((f) => f.path === 'system/fvSolution').text), 'geometric interface variants');
  const p = built.openfoamPipeCase, g = p.files.find((f) => f.path === 'constant/g').text.match(/value\s+\(([-\d.e]+) ([-\d.e]+) ([-\d.e]+)\)/).slice(1).map(Number);
  near(Math.hypot(...g), 9.80665, 1e-6, 'gravity magnitude'); near(-g[0] / 9.80665, Math.sin((snap.angleDeg * Math.PI) / 180), 1e-4, 'gravity resolved along the pipe axis');
  const U = p.files.find((f) => f.path === '0.orig/U').text, q = [...U.matchAll(/volumetricFlowRate ([\d.e-]+);/g)].map((m) => +m[1]), A = (Math.PI / 4) * snap.D ** 2;
  near(q[0], snap.vsl * A, 1e-5, 'liquid inlet flow rate = vsl × area'); near(q[1], snap.vsg * A, 1e-5, 'gas inlet flow rate = vsg × area');
  const tp = p.files.find((f) => f.path === 'constant/transportProperties').text; ok(tp.includes(`rho             ${B.ff(snap.rhoL, 6)}`) && tp.includes(`sigma           ${B.ff(snap.sigma, 5)}`), 'fluid properties of the location are written');
  ok(B.checkBlockMesh(p.mesh).ok && p.plan.cells === p.mesh.nCells && p.plan.coreHours > 0 && p.plan.yPlus > 0, 'pipe case: mesh passes its checks, plan has cells and cost');
  near(B.blockMeshVolume(p.mesh), (Math.PI / 4) * snap.D ** 2 * 20 * snap.D, 3e-4, 'pipe case: meshed volume = π/4 D² L');
  const fine = B.pipeMeshPlan(snap, { turbulence: 'WALE', cellsPerDiameter: 60, yPlus: 1 }), coarse = B.pipeMeshPlan(snap, { turbulence: 'kOmegaSST', cellsPerDiameter: 20 });
  ok(fine.cells > 10 * coarse.cells && fine.coreHours > coarse.coreHours && fine.growth <= 1.2501 && fine.dns.cells > fine.cells, 'mesh plan: finer target costs more; DNS estimate exceeds the LES mesh');
  ok(/DNS/.test(B.openfoamPipeCase({ case: snap, turbulence: 'DNS' }).readme) && B.pipeMeshPlan(snap, { turbulence: 'DNS' }).notes.some((n) => /core-hours/.test(n)), 'DNS case warns about the required resolution and cost');
  const bend = B.openfoamPipeCase({ case: snap, geometry: 'bend', cellsPerDiameter: 12 }); ok(B.checkBlockMesh(bend.mesh).ok && bend.mesh.stations.some((s) => s.kind === 'bend') && bend.reference.forceN > 0, 'pipe case with a 90° bend');
  for (const [dim, solver] of [['2d', 'compressibleInterFoam'], ['3d', 'interFoam']]) { const r = B.openfoamRiserSlugCase({ case: snap, dimension: dim, solver }); ok(B.checkBlockMesh(r.mesh).ok && r.mesh.type === (dim === '2d' ? 'slot' : 'ogrid'), `riser case ${dim}: mesh`); const e = r.mesh.stations[r.mesh.stations.length - 1].t, gg = r.files.find((f) => f.path === 'constant/g').text.match(/value\s+\(([-\d.e]+) ([-\d.e]+) ([-\d.e]+)\)/).slice(1).map(Number); near((e[0] * gg[0] + e[1] * gg[1] + e[2] * gg[2]) / 9.80665, -1, 1e-6, `riser case ${dim}: the riser points straight up against gravity`); }
  const cp = built.openfoamCoupledCase; ok(cp.synthetic && cp.series.t.length > 8 && /volumetricFlowRate table/.test(cp.files.find((f) => f.path === '0.orig/U').text), 'coupled case: slug-train inlet tables without a transient result');
  { const s = cp.series, T = s.t[s.t.length - 1] - s.t[0]; let iL = 0; for (let i = 1; i < s.t.length; i++) iL += 0.5 * (s.vsl[i] + s.vsl[i - 1]) * (s.t[i] - s.t[i - 1]); near(iL / T, snap.vsl, 0.08, 'coupled case: slug train keeps the mean liquid rate'); }
  const withSeries = B.openfoamCoupledCase({ case: { ...snap, series: { t: [0, 10, 20, 30, 40], vsl: [1, 2, 1, 2, 1], vsg: [2, 1, 2, 1, 2] } } }); ok(!withSeries.synthetic && withSeries.plan.endTime === 40, 'coupled case: uses the 1-D transient series when the case has one');
  for (const geometry of ['span', 'bend', 'jumper']) for (const element of ['shell', 'solid']) for (const ends of ['fixed', 'pinned']) { const r = B.calculixPipeCase({ case: snap, geometry, element, ends, nC: 8, elemLen: 0.6 }), c = B.checkFeMesh(r.mesh), st = r.files.find((f) => f.path === 'pipe_static.inp').text, dy = r.files.find((f) => f.path === 'pipe_dynamic.inp').text; ok(c.ok && /\*STATIC/.test(st) && /\*TEMPERATURE/.test(st) && /GRAV/.test(st) && /\*FREQUENCY/.test(r.files.find((f) => f.path === 'pipe_modal.inp').text) && /\*AMPLITUDE, NAME=SLUG/.test(dy) && /\*CLOAD, AMPLITUDE=SLUG/.test(dy) && (element === 'shell' ? /\*DYNAMIC, DIRECT/ : /\*MODAL DYNAMIC/).test(dy), `CalculiX ${geometry} ${element} ${ends}: static, thermal, modal and dynamic steps`); }
  { const r = built.calculixPipeCase, ref = r.reference, ri = snap.D / 2, ro = ri + snap.wt; near(ref.hoop, ((snap.P * ri - snap.pExt * ro) * 1e5) / snap.wt, 1e-9, 'CalculiX: thin-wall hoop stress hand value'); ok(r.files.find((f) => f.path === 'pipe_static.inp').text.includes(`Eall, P, ${B.ff(((snap.P * ri - snap.pExt * ro) * 1e5) / (ri + snap.wt / 2), 8)}`), 'CalculiX: net pressure on the shell mid-surface'); }
  { const dy = B.calculixPipeCase({ case: snap, geometry: 'bend', element: 'solid' }), t = dy.files.find((f) => f.path === 'pipe_dynamic.inp').text, load = [...t.matchAll(/^BEND, (\d), ([-\d.e]+)$/gm)].map((m) => +m[2]), n = dy.mesh.sets.BEND.length; near(Math.hypot(...load) * n, dy.reference.forceN, 1e-6, 'CalculiX: nodal slug loads add up to the peak slug force'); near(dy.reference.forceN, snap.slug.forceN, 1e-9, 'CalculiX: 90° bend carries the √2 ρ A v² slug force of the case'); }
  const fsi = built.fsiCase, xml = fsi.files.find((f) => f.path === 'precice-config.xml').text;
  ok((xml.match(/<[a-zA-Z][^>]*[^/]>/g) || []).length === (xml.match(/<\/[^>]+>/g) || []).length + 1 /* the xml declaration */ && /serial-implicit/.test(xml) && /serial-explicit/.test(B.fsiCase({ case: snap, coupling: 'explicit' }).files[0].text), 'preCICE configuration: tags balance; implicit and explicit schemes');
  ok(B.checkBlockMesh(fsi.fluid).ok && B.checkFeMesh(fsi.solid).ok && fsi.solid.type === 'C3D8I' && fsi.files.find((f) => f.path === 'solid-calculix/mesh.inc').text.includes('NSET=Ninterface'), 'FSI: fluid and solid meshes of the same path, interface node set');
  const py = built.coolpropScript.files[0].text, py2 = built.teqpScript.files[0].text;
  ok(/AbstractState\("HEOS"/.test(py) && py.includes('"Methane"') && py.includes('n-Decane') && (py.match(/^    \(\d+, \d+, /gm) || []).length === 6, 'CoolProp script: HEOS mixture, components mapped, one row per grid point');
  ok(/"kind": "PCSAFT"/.test(py2) && (py2.match(/"sigma_Angstrom"/g) || []).length === grid.comps.length + 1 && /"kind": "CPA"/.test(B.teqpScript({ grid, model: 'CPA' }).files[0].text), 'teqp script: PC-SAFT coefficients for every component; CPA variant');
  near(B.pcsaftPseudo(142.285)[0], 4.6627, 0.04, 'PC-SAFT pseudo-component correlation reproduces n-decane (segment number)'); near(B.pcsaftPseudo(142.285)[2], 243.87, 0.02, 'PC-SAFT pseudo-component correlation reproduces n-decane (energy)');
}

// ---------------------------------------------------------------------------------------------- importers
section('importers on solver output');
{ // interFoam pipe section (real run of the generated case)
  const tabs = B.parseFoamTables(SAMPLE.pipe);
  ok(tabs.length === 5 && tabs.map((t) => t.object).join() === 'forces,holdupVolume,probes,section1,section3', 'collected file is split into its function-object tables');
  ok(tabs[0].columns.join() === 'total_x,total_y,total_z,pressure_x,pressure_y,pressure_z,viscous_x,viscous_y,viscous_z' && tabs[2].columns.join() === 'p.1,p.2,p.3' && tabs[3].columns.join() === 'areaAverage(alpha.liquid),areaAverage(p)', 'column names from the headers (forces, probes, surfaceFieldValue)');
  near(tabs[0].data[0][0], 5.1994189e-1, 1e-9, 'first force value'); near(tabs[3].t[1], 0.0025, 1e-12, 'time column');
  const L = 6 * D, r = B.importOpenfoamPostProcessing(SAMPLE.pipe, { sectionDistance: 0.5 * L, wallArea: Math.PI * D * L }, { discard: 0 }), m = r.metrics;
  near(m.holdupMean, 0.55125, 1e-4, 'mean holdup of the volume average'); ok(m.slugFrequency === 0, 'no slugs in a few milliseconds of start-up');
  ok(m.pressureGradient > 100 && m.pressureGradient < 2000 && m.wallShear > 0.1 && m.wallShear < 20 && m.forcePeak > 300 && m.forcePeak < 600, `pressure gradient ${m.pressureGradient?.toFixed(0)} Pa/m, wall shear ${m.wallShear?.toFixed(2)} Pa and force ${m.forcePeak?.toFixed(0)} N are of the right order`);
  ok(r.series.some((s) => s.name === 'forces: |total|') && r.t.length === tabs[1].t.length, 'series and primary time base');
  const cmp = B.compareFlow(m, { holdup: 0.5755, dpdx: 400, tauW: 5.9 }), h = cmp.find((c) => /volume average/.test(c.quantity));
  near(h.ratio, m.holdupMean / 0.5755, 1e-12, 'comparison ratio'); near(h.difference, m.holdupMean - 0.5755, 1e-9, 'comparison difference'); ok(cmp.every((c) => c.external !== null) && cmp.some((c) => c.unit === 'Pa/m'), 'comparison rows only for imported quantities');
  // the same tables as separate files
  const parts = SAMPLE.pipe.split(/^### file: (.*)$/m).slice(1), sep = []; for (let i = 0; i < parts.length; i += 2) sep.push({ name: parts[i], text: parts[i + 1] });
  near(B.importOpenfoamPostProcessing(sep, { sectionDistance: 0.5 * L, wallArea: Math.PI * D * L }, { discard: 0 }).metrics.holdupMean, m.holdupMean, 1e-12, 'separate files give the same result as the collected file');
}
{ // compressibleInterFoam pipeline–riser (real run, 22 s, decimated): the riser fills and the top section churns
  const r = B.importOpenfoamPostProcessing(SAMPLE.riser, {}), m = r.metrics;
  ok(m.holdupSection > 0.8 && m.holdupSection < 0.9, `riser-base holdup ${m.holdupSection?.toFixed(3)}`); ok(m.holdupMax - m.holdupMin > 0.3, 'riser-top holdup swings');
  ok(m.pressureMean > 8.55e6 && m.pressureMean < 8.65e6 && m.duration > 20, 'riser-base pressure level and duration');
}
{ // DPMFoam erosion case (real run)
  const r = B.importOpenfoamPostProcessing(SAMPLE.erosion, { wallFaceArea: 2e-3 }), m = r.metrics;
  ok(m.erodedVolume > 0 && m.erodedVolumeMax > 0 && m.erodedVolume >= m.erodedVolumeMax, 'eroded volume: total and worst face');
  ok(m.erosionRateMmY > 0 && Number.isFinite(m.erosionRateMmY), `erosion rate ${m.erosionRateMmY?.toExponential(2)} mm/y from the slope of the worst-face series`);
  near(B.importOpenfoamPostProcessing(SAMPLE.erosion, { wallFaceArea: 1e-3 }).metrics.erosionRateMmY, 2 * m.erosionRateMmY, 1e-9, 'erosion rate scales with 1 / face area');
  ok(m.pressureDrop > 1000 && m.pressureDrop < 10000 && r.series.some((s) => s.name === 'cloudInfo: nParcels'), 'pressure drop and parcel count series');
}
{ // slug frequency of a known signal: 0.2 Hz square-ish holdup at the outlet, start-up excluded
  const t = [], y = []; for (let i = 0; i <= 2000; i++) { const ti = i * 0.05; t.push(ti); y.push(0.3 + 0.5 * (Math.sin(2 * Math.PI * 0.2 * ti) > 0.6 ? 1 : 0) + 0.01 * Math.sin(37 * ti)); }
  const f = B.fluctuation(t, y); near(f.freqCrossing, 0.2, 0.02, 'slug frequency by threshold crossings'); near(f.freqSpectrum, 0.2, 0.06, 'slug frequency by spectrum peak');
  const csv = 'time,alpha.liquid\n' + t.map((ti, i) => `${ti},${y[i]}`).join('\n'), r = B.importOpenfoamPostProcessing([{ name: 'postProcessing/outletHoldup/0/series.csv', text: csv }], {});
  near(r.metrics.slugFrequency, 0.2, 0.02, 'slug frequency through the importer (CSV with a header)');
  ok(B.fluctuation([0, 1, 2, 3, 4, 5, 6, 7, 8, 9], new Array(10).fill(0.4)).freqCrossing === 0, 'a constant signal has no frequency');
}
{ // CalculiX .dat and .frd (real runs of the generated decks)
  const md = B.importCalculixDat(SAMPLE.modalDat);
  ok(md.frequencies.length === 10, 'eigenvalue table: ten modes'); near(md.metrics.f1, 34.60509, 1e-6, 'first natural frequency (Hz)'); near(md.frequencies[2], 91.97187, 1e-6, 'third natural frequency');
  near(md.frequencies[2] / md.metrics.f1, (7.8532046 / 4.730040745) ** 2, 0.05, 'second bending frequency over the first is close to the clamped-beam ratio (shell: shear deformation lowers the higher mode)');
  const sd = B.importCalculixDat(SAMPLE.staticDat);
  ok(sd.displacements.length === 1 && sd.displacements[0].set === 'BEND' && sd.displacements[0].n >= 4 && sd.metrics.maxDisplacement > 0, 'displacement block of a node set'); near(sd.forces[0].f[0], 4.124111e5, 1e-9, 'total reaction force of the first step'); ok(sd.metrics.peakReaction >= Math.hypot(...sd.forces[0].f), 'peak reaction over the steps');
  const fr = B.importCalculixFrd(SAMPLE.frd);
  ok(fr.nodes === 4 && fr.steps.length === 4 && fr.steps.map((s) => s.kind).join() === 'DISP,STRESS,DISP,STRESS' && fr.steps[2].time === 2, '.frd: nodes and displacement / stress blocks of both steps');
  ok(fr.metrics.maxMises > 1e6 && fr.metrics.maxMises < 5e8 && fr.metrics.maxDisplacement > 0, '.frd: von Mises stress and displacement maxima');
  { const l = SAMPLE.frd.split('\n').filter((x) => x.startsWith(' -1')), i0 = SAMPLE.frd.split('\n').findIndex((x) => /^ -4  STRESS/.test(x)), row = SAMPLE.frd.split('\n').slice(i0).find((x) => x.startsWith(' -1')), v = [0, 1, 2, 3, 4, 5].map((k) => +row.slice(13 + 12 * k, 25 + 12 * k)), vm = Math.sqrt(0.5 * ((v[0] - v[1]) ** 2 + (v[1] - v[2]) ** 2 + (v[2] - v[0]) ** 2) + 3 * (v[3] ** 2 + v[4] ** 2 + v[5] ** 2)); ok(l.length > 8 && fr.steps[1].max >= vm * (1 - 1e-12), '.frd: fixed-width stress columns (values that touch) are read'); }
  const cs = B.compareStructure({ f1: 34.6, maxMises: 160e6, maxDisplacement: 1.8e-4 }, { f1: 35.64, smys: 448e6 });
  near(cs.find((c) => c.unit === 'Hz').ratio, 34.6 / 35.64, 1e-12, 'structure comparison: frequency ratio'); near(cs.find((c) => /SMYS/.test(c.quantity)).external, 160 / 448, 1e-12, 'structure comparison: utilisation');
}
{ // untrusted input
  let n = 0; for (const bad of [() => B.parseFoamTables('hello world'), () => B.importCalculixFrd('not a result file'), () => B.importOpenfoamPostProcessing(42), () => B.importPropertyTable('a,b\n1,2\n'), () => B.importCalculixDat({}), () => B.parseFoamTables('x'.repeat(B.BRIDGE_LIMITS.chars + 1))]) { try { bad(); } catch (e) { if (e instanceof Error && e.message.length > 5) n++; } }
  ok(n === 6, 'malformed or oversized input is refused with a message');
  const odd = B.parseFoamTables('# Time a b\n0 1 2\n0.1 nan 3\n0.2 (4 5)\n0.3 6 7\n__proto__ 1 2\n'); ok(odd[0].t.join() === '0,0.2,0.3' && odd[0].data[0].join() === '1,4,6', 'non-numeric rows are skipped, parentheses stripped');
  const rec = B.bridgeRecord('openfoamPipeCase', 'OpenFOAM', { a: 1.23456789, b: NaN, c: 'x'.repeat(500), d: {} }, [{ quantity: 'q', unit: 'u', external: 1, inApp: null, ratio: Infinity, difference: undefined }]);
  ok(rec.metrics.a === 1.23457 && !('b' in rec.metrics) && rec.metrics.c.length === 80 && !('d' in rec.metrics) && rec.comparison[0].inApp === null && JSON.stringify(rec).length < 800, 'stored record is small and clean');
}

// ---------------------------------------------------------------------------------------------- property tables
section('property tables');
{
  const table = buildTable(DEFAULT_FLUID, { nP: 5, nT: 4 }), csv = B.propertyTableCSV(table), back = B.importPropertyTable(csv, table);
  ok(csv.split('\n')[0] === 'P,T,' + B.TABLE_FIELDS.join(',') && csv.trim().split('\n').length === 21, 'CSV has one row per grid point');
  ok(back.missing === 0 && back.replaced.length === B.TABLE_FIELDS.length, 'round trip fills every field');
  ok(JSON.stringify(Object.keys(back.table).sort()) === JSON.stringify(Object.keys(table).sort()), 'imported table has exactly the keys of buildTable');
  ok(B.TABLE_FIELDS.every((f) => JSON.stringify(back.table[f]) === JSON.stringify(table[f])) && JSON.stringify(back.table.P) === JSON.stringify(table.P) && JSON.stringify(back.table.lnP) === JSON.stringify(table.lnP) && JSON.stringify(back.table.spec) === JSON.stringify(table.spec) && JSON.stringify(back.table.rates) === JSON.stringify(table.rates), 'round trip is exact');
  near(lookup(back.table, 37, 55).rhoG, lookup(table, 37, 55).rhoG, 1e-14, 'kernel lookup works on the imported table');
  ok(fluidModel({ fluid: DEFAULT_FLUID, outputs: { pvt: { table: back.table } } }).table === back.table, 'the fluid model accepts it as the published PVT table');
  const alone = B.importPropertyTable(csv); ok(alone.table.rhoG[2][1] === table.rhoG[2][1] && alone.table.eosId === 'imported', 'a complete table imports without a base table');
  // partial table (as the scripts write): only gas columns, some cells empty
  const part = ['P,T,rhoG,zG,source', ...table.P.flatMap((p, i) => table.T.map((t, j) => `${p},${t},${i === 0 && j === 0 ? '' : table.rhoG[i][j] * 1.01},${i === 0 && j === 0 ? '' : table.zG[i][j]},model X`))].join('\n'), imp = B.importPropertyTable(part, table);
  ok(imp.replaced.join() === 'rhoG,zG' && imp.table.rhoG[0][0] === table.rhoG[0][0] && imp.source === 'model X', 'partial table: empty cells and absent fields come from the case table');
  near(imp.table.rhoG[3][2], 1.01 * table.rhoG[3][2], 1e-12, 'partial table: supplied cells replace the case values');
  const dev = B.comparePropertyTables(table, imp.table).find((d) => d.field === 'rhoG'); near(dev.meanPct, 100 * (1 - 1 / 1.01), 1e-6, 'deviation statistics'); ok(dev.bias < 0, 'bias sign: in-app below the imported values');
  let t1 = false, t2 = false, t3 = false; try { B.importPropertyTable('P,T,rhoG\n1,10,5\n2,10,6\n1,20,5\n2,20,6\n', table); } catch { t1 = true; } try { B.importPropertyTable(part); } catch { t2 = true; } try { B.importPropertyTable(part.replace(/,model X/g, ',x').replace(String(table.rhoG[2][2] * 1.01), '-5'), table); } catch { t3 = true; }
  ok(t1 && t2 && t3, 'a different grid, an incomplete table without a base, and non-physical values are refused');
  // real GERG-2008 rows written by the generated CoolProp script (CoolProp 6.7.0)
  const g = SAMPLE.gergHead.trim().split('\n'); ok(g[0] === 'P,T,wG,rhoG,zG,cpG,jtG,wSound,rhoO,cpO,source' && /GERG-2008 \(CoolProp/.test(g[1]), 'CoolProp script output has the columns the importer expects');
  const pg = B.propertyGrid(DEFAULT_FLUID, { P: [20, 100], T: [20, 80] }); ok(pg.points.length === 4 && pg.points.every((p) => Math.abs(p.x.reduce((a, b) => a + b, 0) - 1) < 1e-9 && Math.abs(p.y.reduce((a, b) => a + b, 0) - 1) < 1e-9 && p.beta >= 0 && p.beta <= 1) && pg.comps.length === pg.points[0].x.length, 'property grid: normalised phase compositions at every point');
}
if (process.env.BRIDGE_PYTHON) { // execute the generated scripts in the real libraries and re-import their output
  section('property scripts executed with ' + process.env.BRIDGE_PYTHON);
  const fm = fluidModel({}), full = B.propertyGrid(fm.spec, { P: fm.table.P.filter((_, i) => i % 3 === 0), T: fm.table.T.filter((_, j) => j % 4 === 0) }), base = buildTable(fm.spec, { nP: 8, nT: 5 });
  base.P = full.P.slice(); base.T = full.T.slice(); base.lnP = base.P.map(Math.log); for (const f of B.TABLE_FIELDS) base[f] = full.P.map((_, i) => full.T.map((__, j) => fm.table[f][i * 3][j * 4]));
  const dir = fs.mkdtempSync(path.join(os.tmpdir(), 'hsbridge-'));
  for (const [r, lim] of [[B.coolpropScript({ grid: full }), 5], [B.teqpScript({ grid: full }), 8]]) {
    const f = path.join(dir, r.files[0].path); fs.writeFileSync(f, r.files[0].text);
    try { execFileSync(process.env.BRIDGE_PYTHON, [f], { cwd: dir, stdio: ['ignore', 'pipe', 'pipe'], timeout: 600000 }); const imp = B.importPropertyTable(fs.readFileSync(path.join(dir, r.output), 'utf8'), base), dev = B.comparePropertyTables(base, imp.table).find((d) => d.field === 'rhoG'); ok(imp.filled > 20 && dev && dev.meanPct < lim, `${r.name}: runs, re-imports, gas density of the cubic EOS within ${lim} % on average (${dev?.meanPct.toFixed(2)} %)`); }
    catch (e) { ok(false, `${r.name}: ${String(e.stderr || e.message).slice(0, 300)}`); }
  }
  fs.rmSync(dir, { recursive: true, force: true });
}

// ---------------------------------------------------------------------------------------------- page
section('page (DOM stub)');
{
  class Node { constructor() { this.parent = null; } remove() { const p = this.parent; if (p) p.children = p.children.filter((c) => c !== this); this.parent = null; } }
  class Text extends Node { constructor(t) { super(); this.nodeType = 3; this.textContent = String(t); } walk() {} }
  class El extends Node {
    constructor(t) { super(); this.tagName = String(t).toUpperCase(); this.children = []; this.attrs = {}; this.dataset = {}; this.listeners = {}; this.nodeType = 1; this.className = ''; this.classList = { add() {}, remove() {}, toggle() {}, contains: () => false }; }
    append(...a) { for (const x of a) { x.parent = this; this.children.push(x); } } get firstChild() { return this.children[0] || null; }
    setAttribute(k, v) { this.attrs[k] = String(v); } getAttribute(k) { return this.attrs[k] ?? null; } addEventListener(ev, fn) { (this.listeners[ev] = this.listeners[ev] || []).push(fn); }
    click() { for (const fn of this.listeners.click || []) fn({ target: this }); } scrollIntoView() {} focus() {} querySelector() { return null; } toBlob() {} getBoundingClientRect() { return { left: 0, top: 0, width: 300, height: 200 }; }
    get isConnected() { return false; } get clientWidth() { return 0; }
    walk(fn) { fn(this); for (const c of this.children) c.walk(fn); }
    get textContent() { return this._text ?? this.children.map((c) => c.textContent).join(''); } set textContent(v) { this._text = String(v); this.children = []; }
  }
  globalThis.document = { createElement: (t) => new El(t), createTextNode: (t) => new Text(t), body: new El('body') };
  globalThis.ResizeObserver = class { observe() {} disconnect() {} };
  const saved = []; const origCreate = URL.createObjectURL; URL.createObjectURL = (b) => { saved.push(b); return 'blob:x'; }; URL.revokeObjectURL = () => {};
  const { bridgePage, BRIDGE_CARDS } = await import('../js/pages/bridge.js');
  const root = document.createElement('div');
  let err = null; try { bridgePage(root, {}); } catch (e) { err = e; }
  ok(!err, 'page renders with nothing solved yet' + (err ? ': ' + (err.stack || err.message).split('\n').slice(0, 3).join(' | ') : ''));
  const all = []; root.walk((e) => all.push(e));
  const cards = all.filter((e) => e.tagName === 'SECTION' && e.dataset.generator), buttons = all.filter((e) => e.tagName === 'BUTTON'), text = root.textContent;
  ok(cards.length === BRIDGE_CARDS.length && cards.length >= 10, `one card per hand-off (${cards.length})`);
  ok(BRIDGE_CARDS.every((c) => typeof B[c.generator] === 'function') && Object.values(B.HANDOFF).flat().every((e) => BRIDGE_CARDS.some((c) => c.generator === e.generator)), 'every hand-off entry has a card on the page');
  ok(/External solvers/.test(text) && /ready-to-run case/.test(text) && /Commands to run/.test(text) && /\.\/Allrun/.test(text) && /Estimated cost/.test(text), 'explanation, commands and cost estimates are shown');
  ok(all.filter((e) => e.tagName === 'INPUT' && e.attrs.type === 'range').length === 1 && all.filter((e) => e.tagName === 'INPUT' && e.attrs.type === 'file' && 'multiple' in e.attrs).length === cards.length, 'location slider and one import control per card');
  ok(all.every((e) => !('style' in e.attrs)) && all.filter((e) => e.tagName === 'INPUT' && !e.hidden && e.attrs.type !== 'file').every((e) => e.attrs.id || e.attrs['aria-label']), 'no inline styles; inputs are labelled');
  const dls = buttons.filter((b) => /Download case/.test(b.textContent)); let dlErr = null; try { dls.forEach((b) => b.click()); } catch (e) { dlErr = e; }
  ok(!dlErr && dls.length === cards.length && saved.length === cards.length, 'every card writes its .zip' + (dlErr ? ': ' + dlErr.message : ''));
  if (saved.length) { const z = B.unzipStored(new Uint8Array(await saved[0].arrayBuffer())); ok(z.length > 15 && z.every((e) => e.crcOk) && z.some((e) => /\/system\/blockMeshDict$/.test(e.path)), 'downloaded archive is a complete case'); }
  URL.createObjectURL = origCreate;
}

console.log(`\n${count - fails} of ${count} checks passed.`);
if (fails) { console.log(`${fails} FAILED`); process.exit(1); }
