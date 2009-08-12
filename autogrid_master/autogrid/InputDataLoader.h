/*
    FastGrid (formerly AutoGrid)

    Copyright (C) 2009 The Scripps Research Institute. All rights reserved.
    Copyright (C) 2009 Masaryk University. All rights reserved.

    AutoGrid is a Trade Mark of The Scripps Research Institute.

    This program is free software; you can redistribute it and/or
    modify it under the terms of the GNU General Public License
    as published by the Free Software Foundation; either version 2
    of the License, or (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program; if not, write to the Free Software
    Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.
*/

#pragma once
#include "GridMap.h"

class InputDataLoader : public InputData
{
public:
    InputDataLoader(LogFile *logFile);
    ~InputDataLoader();
    void load(const char *gridParameterFilename, GridMapList &gridmaps, ParameterLibrary &parameterLibrary);

private:
    LogFile *logFile;

    int checkSize(int numGridPointsMinusOne, char axischar);
    static int parseGPFLine(const char *line);
    static int parseTypes(char *line, char **words, int maxwords);
    static int strIndex(const char *s, const char *t);
};
