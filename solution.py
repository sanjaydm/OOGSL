x = [
[0,	0],
[0,	0.000173444947667069],
[3.58216705473186e-05,	3.51685502260349e-05],
[0.000329173276856392,	0.000155698586523699],
[0.000138060535373003,	6.64165992713399e-05],
[0.000622644196069521,	0.00013774267088525],
[0.000298872111139876,	9.37407291312628e-05],
[0.000880023849192414,	0.000119611158198887],
[0.000510405451045038,	0.000117137877278268],
[0.00110098452939508,	0.000101329846404726],
[0.000764803793014322,	0.000136605295054052],
[0.00128525883795634,	8.29297792660518e-05],
[0.00105420515318171,	0.000152140503128192],
[0.00143264086515164,	6.44405486181741e-05],
[0.0013707428230394,	0.000163741246305651],
[0.00154297223649449,	4.58830189755961e-05],
[0.00170654581849576,	0.000171405514657442],
[0.00161614675689252,	2.728675940573e-05],
[0.00205373934735143,	0.000175131494119282],
[0.00165210667026489,	8.67497824216424e-06],
[0.00240444518995367,	0.000174917569813819],
[0.00165084290849132,	-9.93472958432666e-06],
[0.00275078205853332,	0.000170762301445421],
[0.00161238354168842,	-2.85179019390722e-05],
[0.00308486602216148,	0.000162664477677248],
[0.00153680100420332,	-4.70562482475916e-05],
[0.00339881088453083,	0.000150623043603983],
[0.0014241975043697,	-6.55331655889126e-05],
[0.00368472859170436,	0.00013463720444801],
[0.00127471968991543,	-8.39309633509568e-05],
[0.00393472969177223,	0.000114706357257935],
[0.00108853812452611,	-0.000102234881403318],
[0.00414092378257738,	9.08301589892452e-05],
[0.000865857694454108,	-0.00012042894642712],
[0.00429541999127685,	6.30084816238848e-05],
[0.000606900678981775,	-0.00013850374095688],
[0.00439032754382852,	3.12415579088446e-05],
[0.000311924713353213,	-0.000156449633731736],
[0.0044177563837335,	-4.47011196536499e-06],
[-1.88030853371089e-05,	-0.000174253017141541],
[0.00436981786503532,	-4.41256221842937e-05],
[-0.00038498796620386,	-0.000191908771698864],
[0.00423862550958385,	-8.77237043829768e-05],
[-0.000786329678898472,	-0.000209412898322007],
[0.00401629590547502,	-0.000135262578857178],
[-0.0012225050053875,	-0.000226780597514572],
[0.00369494974174667,	-0.000186740145838208],
[-0.00169322053949956,	-0.000244166372240338],
[0.00326671289107632,	-0.000242154658330005],
[-0.00219836130849547,	-0.000262534219658263],
[0.00272371762919541,	-0.000301509785291212],
[-0.00273901636060059,	-0.000287848851064675],
[0.00206349040406362,	-0.00035398529674288],
[-0.00320918453966376,	-0.000158946803103605],
[0.00134281645817681,	-0.000356555603251353],
[-0.00323394242677737,	0.000124541970760594],
[0.000677734757040755,	-0.000298405563018702],
[-0.00270649776974218,	0.0004013493571451],
[0.000189663129144342,	-0.0001795490770392],
[-0.00162799509071739,	0.000676796593896415],
[0,	0],
[0,	0.000950876336601207],]
from matplotlib import pyplot as plt
from numpy import *
x = matrix(x)
plt.plot(x[0:len(x):2],'*-')
plt.show()
