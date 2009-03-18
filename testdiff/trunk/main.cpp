#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <sstream>
#include <cmath>

#define implication(x, y) (!(x) || (y))

// Max
template<typename T>
inline T max(T x, T y)
{
    return x > y ? x : y;
}

template<typename T>
inline std::string tostr(const T &var)
{
    std::ostringstream stream;
    stream << var;
    return stream.str();
}

template<typename T>
inline bool strto(const std::string &str, T &out)
{
    std::istringstream stream(str);
    stream >> out;
    return !!stream;
}

class FileDiff
{
public:
    FileDiff()
    {
        fileinfo.reserve(1024);
        added.reserve(256);
        removed.reserve(256);
    }

    void ReadFileInfo(const std::string &command, std::istream &input);
    void InsertLine(const std::string line);
    void Clear() { added.clear(); removed.clear(); }
    bool ExamineChanges(std::ostream &output, double &absErrorMax, double &relErrorMax);

private:
    std::string line, fileinfo;
    std::vector<std::string> added;
    std::vector<std::string> removed;

    static bool IsSimilar(const std::string &s1, const std::string &s2, double &absErrorMax, double &relErrorMax);

    friend std::ostream &operator <<(std::ostream &s, const FileDiff &diff);
};

std::ostream &operator <<(std::ostream &s, const FileDiff &diff)
{
    s << diff.fileinfo << '\n';
    int i = 1;
    for (std::vector<std::string>::const_iterator it = diff.removed.begin(); it != diff.removed.end(); ++it, ++i)
        s << '[' << i << "]-" << *it << '\n';
    i = 1;
    for (std::vector<std::string>::const_iterator it = diff.added.begin(); it != diff.added.end(); ++it, ++i)
        s << '[' << i << "]+" << *it << '\n';
    return s;
}

void FileDiff::ReadFileInfo(const std::string &command, std::istream &input)
{
    char buf[512];
    fileinfo = command;
    fileinfo += "\n";
    input.getline(buf, 512);
    fileinfo += buf;
    fileinfo += "\n";
    input.getline(buf, 512);
    fileinfo += buf;
}

void FileDiff::InsertLine(const std::string line)
{
    if (line[0] == '+')
        added.push_back(std::string(line, 1));
    else if (line[0] == '-')
        removed.push_back(std::string(line, 1));
}

// examines changes and returns true if they were irrelevant (typically caused by a different floating-point implementation)
bool FileDiff::ExamineChanges(std::ostream &output, double &absErrorMax, double &relErrorMax)
{
    if (!added.size() && !removed.size())
        return true;

    if (added.size() != removed.size())
    {
        output << "*** Test failed. The count of '+'s is different from the count of '-'s. See below. ***\n" << *this;
        return false;
    }

    bool ok = true;
    for (size_t i = 0; i < added.size(); i++)
    {
        double absError = 0, relError = 0;
        if (!IsSimilar(added[i], removed[i], absError, relError))
        {
            if (ok)
            {
                output << "*** Test failed. See below. ***\n";
                output << fileinfo << '\n';
            }
            output << "        [" << (i+1) << "] -{" << removed[i] << "} +{" << added[i] << "}, AE: " << absError << ", RE: " << relError << "\n";
            ok = false;
        }

        absErrorMax = max(absErrorMax, absError);
        relErrorMax = max(relErrorMax, relError);
    }
    return ok;
}

bool FileDiff::IsSimilar(const std::string &s1, const std::string &s2, double &absErrorMax, double &relErrorMax)
{
    // do not care about comments
    if (*s1.c_str() == '#' && *s2.c_str() == '#')
        return true;

    // do not care about execution times
    if (s1.find("Real") != std::string::npos && s2.find("Real") != std::string::npos &&
        s1.find("CPU") != std::string::npos && s2.find("CPU") != std::string::npos &&
        s1.find("System") != std::string::npos && s2.find("System") != std::string::npos)
        return true;

    // do not care about path differences
    if (s1.find("../") != std::string::npos || s2.find("../") != std::string::npos)
        return true;

    // do not care about the date and hostname
    if (s1.find("This file was created at:") != std::string::npos && s2.find("This file was created at:") != std::string::npos)
        return true;
    if (s1.find("using:") != std::string::npos && s2.find("using:") != std::string::npos)
        return true;

    std::istringstream stream1(s1), stream2(s2);
    std::string n1, n2;
    double x,y;

    // both strings may consist of several separate numbers or substrings, we must check each of them
    while (!stream1.eof() || !stream2.eof())    // until both of them reach eof
    {
        stream1 >> n1;
        stream2 >> n2;

        double absError = 0;
        double relError = 0;

        // a parser error occured
        if (!stream1 || !stream2)
            return false;

        // substrings are the same
        if (n1 == n2)
            continue;

        // convert both strings to double
        bool converted = strto<double>(n1, x) && strto<double>(n2, y);
        if (!converted)
            return false;

        if (x != y)
        {
            absError = fabs(x - y);
            relError = fabs((x - y) / y);
        }

        absErrorMax = max(absErrorMax, absError);
        relErrorMax = max(relErrorMax, relError);

        if (!(0
                // take the absolute error into account
                || (absError < 0.05)
                // take the relative error into account
                || (relError < 0.08)
            ))
            return false;
    }
    return true;
}

///////////////////////////////////////////////////////////////////////////////

int main(int argc, char **argv)
{
    if (argc != 3)
    {
        std::cout << "usage: testdiff [diff_filename] [autogrid_log]\n";
        return 0;
    }

    std::ifstream f;
    f.open(argv[1]);
    if (!f.is_open())
    {
        std::cout << "cannot open file: " << argv[1] << "\n";
        return 1;
    }

    int returnValue = 0;
    char buf[512];
    std::string line;
    line.reserve(256);
    FileDiff diff;
    double absErrorMax = 0, relErrorMax = 0;

    bool skipLog = false;
    while (f.getline(buf, 512))
    {
        line = buf;

        switch (line[0])
        {
        case 'd':
            if (!diff.ExamineChanges(std::cout, absErrorMax, relErrorMax))
                returnValue = 1;
            skipLog = line.find(argv[2]) != std::string::npos;
            if (!skipLog)
            {
                diff.Clear();
                if (std::string(line, 0, 4) == "diff")
                    diff.ReadFileInfo(line, f);
                else
                {
                    std::cerr << "Error: unknown identifier '" << std::string(line, 0, line.find(' ')) << "'.\n";
                    return 2;
                }
            }
            break;

        case '+':
        case '-':
            if (!skipLog)
                diff.InsertLine(line);
            break;
        }
    }

    if (!diff.ExamineChanges(std::cout, absErrorMax, relErrorMax))
        returnValue = 1;

    std::cout << "*** Max AbsError: " << absErrorMax << ", Max RelError: " << relErrorMax << " ***\n";
    if (!returnValue)
        std::cout << "*** Test OK. ***\n";
    return returnValue;
}
