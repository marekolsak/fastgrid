#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <sstream>
#include <cmath>

#define implication(x, y) (!(x) || (y))

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
    bool ExamineChanges(std::ostream &output);

private:
    std::string line, fileinfo;
    std::vector<std::string> added;
    std::vector<std::string> removed;

    static bool IsSimilar(const std::string &s1, const std::string &s2);

    friend std::ostream &operator <<(std::ostream &s, const FileDiff &diff);
};

std::ostream &operator <<(std::ostream &s, const FileDiff &diff)
{
    s << diff.fileinfo << '\n';
    int i = 0;
    for (std::vector<std::string>::const_iterator it = diff.removed.begin(); it != diff.removed.end(); ++it, ++i)
        s << '[' << i << "]-" << *it << '\n';
    i = 0;
    for (std::vector<std::string>::const_iterator it = diff.added.begin(); it != diff.added.end(); ++it, ++i)
        s << '[' << i << "]-" << *it << '\n';
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
bool FileDiff::ExamineChanges(std::ostream &output)
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
        if (!IsSimilar(added[i], removed[i]))
        {
            output << "*** Test failed. The " << (i+1) << (i%10 == 0 && i != 10? "st" : i%10 == 1 && i != 11? "nd" : "th") << " item is different. See below. ***\n";
            ok = false;
        }
    if (!ok)
        output << *this;
    return ok;
}

bool FileDiff::IsSimilar(const std::string &s1, const std::string &s2)
{
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
    while (!stream1.eof() || !stream2.eof())    // until both of them reach an eof
    {
        stream1 >> n1;
        stream2 >> n2;

        // a parser error occured
        if (!stream1 || !stream2)
            return false;

        // substrings are the same
        if (n1 == n2)
            continue;

        // convert both strings to double
        bool converted = strto<double>(n1, x) && strto<double>(n2, y);

        if (!(
                // consider the last digit of the number a rounding error
                (n1.length() == n2.length() && std::string(n1, 0, n1.length()-1) == std::string(n2, 0, n2.length()-1)) ||
                // just an implication
                implication(
                    converted, // implies:
                    // the original string representation may differ even though the numbers are the same
                    (x == y) ||
                    // consider +-0 a rounding error
                    ((x == 0 || y == 0) && (fabs(x) == fabs(y)))
                )
            ))
            return false;
    }
    return true;
}

///////////////////////////////////////////////////////////////////////////////

int main(int argc, char **argv)
{
    if (argc != 2)
    {
        std::cout << "usage: testdiff [filename]\n";
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

    while (f.getline(buf, 512))
    {
        line = buf;

        switch (line[0])
        {
        case 'd':
            if (!diff.ExamineChanges(std::cout))
                returnValue = 1;
            diff.Clear();
            if (std::string(line, 0, 4) == "diff")
                diff.ReadFileInfo(line, f);
            break;

        case '+':
        case '-':
            diff.InsertLine(line);
            break;
        }
    }

    if (!diff.ExamineChanges(std::cout))
        returnValue = 1;

    if (!returnValue)
        std::cout << "*** Test OK. ***\n";
    return returnValue;
}
