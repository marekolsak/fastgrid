#pragma once

// Precalculates the exponential function for receptor and ligand desolvation
class DesolvExpFunc
{
public:
    DesolvExpFunc(double coeffDesolv);
    ~DesolvExpFunc();

    double operator ()(int i) const { return func[i]; }

private:
    double *func;
};
