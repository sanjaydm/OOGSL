x = [0,0.025,0.05,0.075,0.1,0.125,0.15,0.175,0.2,0.225,0.25,0.275,0.3,0.325,0.35,0.375,0.4,0.425,0.45,0.475,0.5,0.525,0.55,0.575,0.6,0.625,0.65,0.675,0.7,0.725,0.75,0.775,0.8,0.825,0.85,0.875,0.9,0.925,0.95,0.975,1,]
y = [-0.0980000698571783,-0.0979999892771813,-0.0980000028082163,-0.097999997969457,-0.0980000022364931,-0.0979999964532148,-0.09800000666253,-0.0979999868503896,-0.0980000254508198,-0.0979999516544763,-0.0980000915191397,-0.0979998267376016,-0.0980003278335092,-0.0979993799671322,-0.0980011727643569,-0.0979977815400283,-0.0980041960081011,-0.0979920648683848,-0.0980150052143547,-0.0979716267909377,-0.0980536492309876,-0.0978985576475364,-0.0981918131699712,-0.0976373065130822,-0.0986858071968146,-0.096703225416269,-0.100452035839868,-0.0933635126897216,-0.106767007748078,-0.0814226610444123,-0.129345865696091,-0.0387266684153551,-0.210089941721336,0.114049869506419,-0.499504872354503,0.665705808421885,-1.56736419031504,2.88045003553585,-6.89799612616703,21.1988697531689,-176.071577358479,]
from matplotlib import pyplot as plt
from numpy import *
x = matrix(x)
y = matrix(y)
xmin = -matrix(x)
y = matrix(y)
plt.plot(x.T, y.T,'*-')
plt.axis('tight')
plt.show()
