#include "desolvexpfunc.h"
#include "autogrid.h"
#include <cmath>

DesolvExpFunc::DesolvExpFunc(double coeffDesolv)
{
    func = new double[MAX_DIST];
    double minusInvTwoSigmaSquared;
    double sigma;

    // exponential function for receptor and ligand desolvation
    // note: the solvation term will not be smoothed
    sigma = 3.6;
    minusInvTwoSigmaSquared = -1 / (2 * sigma * sigma);
    for (int indexR = 1; indexR < MAX_DIST; indexR++)
    {
        double r = angstrom(indexR);
        func[indexR] = exp(sq(r) * minusInvTwoSigmaSquared);
        func[indexR] *= coeffDesolv;
    }
}

DesolvExpFunc::~DesolvExpFunc()
{
    delete [] func;
}
