//Copyright 2015 Ryan Wick

//This file is part of Bandage

//Bandage is free software: you can redistribute it and/or modify
//it under the terms of the GNU General Public License as published by
//the Free Software Foundation, either version 3 of the License, or
//(at your option) any later version.

//Bandage is distributed in the hope that it will be useful,
//but WITHOUT ANY WARRANTY; without even the implied warranty of
//MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//GNU General Public License for more details.

//You should have received a copy of the GNU General Public License
//along with Bandage.  If not, see <http://www.gnu.org/licenses/>.


#include "distance.h"
#include "commoncommandlinefunctions.h"
#include "../program/settings.h"
#include "../graph/assemblygraph.h"
#include "../blast/blastsearch.h"
#include "../program/memory.h"

int bandageDistance(QStringList arguments)
{
    QTextStream out(stdout);
    QTextStream err(stdout);

    if (checkForHelp(arguments))
    {
        printDistanceUsage(&out, false);
        return 0;
    }

    if (checkForHelpAll(arguments))
    {
        printDistanceUsage(&out, true);
        return 0;
    }


    if (arguments.size() < 2)
    {
        printDistanceUsage(&err, false);
        return 1;
    }

    QString graphFilename = arguments.at(0);
    arguments.pop_front();
    if (!checkIfFileExists(graphFilename))
    {
        err << "Bandage error: " << graphFilename << " does not exist." << endl;
        return 1;
    }

    QString queriesFilename = arguments.at(0);
    arguments.pop_front();
    if (!checkIfFileExists(queriesFilename))
    {
        err << "Bandage error: " << queriesFilename << " does not exist." << endl;
        return 1;
    }
    g_settings->blastQueryFilename = queriesFilename;

    //Ensure that the --query option isn't used, as that would overwrite the
    //queries file that is a positional argument.
    if (isOptionPresent("--query", &arguments))
    {
        err << "Bandage error: the --query option cannot be used with Bandage distance." << endl;
        return 1;
    }

    QString error = checkForInvalidDistanceOptions(arguments);
    if (error.length() > 0)
    {
        err << "Bandage error: " << error << endl;
        return 1;
    }

    bool allQueryPaths = false;
    bool sequences = false;
    parseDistanceOptions(arguments, &allQueryPaths, &sequences);

    bool loadSuccess = g_assemblyGraph->loadGraphFromFile(graphFilename);
    if (!loadSuccess)
        return 1;

    if (!createBlastTempDirectory())
    {
        err << "Error creating temporary directory for BLAST files" << endl;
        return 1;
    }

    QString blastError = g_blastSearch->doAutoBlastSearch();
    if (blastError != "")
    {
        err << blastError << endl;
        return 1;
    }

    //Check that there are at least two queries.
    int queryCount = g_blastSearch->m_blastQueries.getQueryCount();
    if (queryCount < 2)
    {
        err << "Bandage error: the queries file must contain at least two sequences." << endl;
        return 1;
    }

    //Report which of the queries were not found.
    QList<BlastQuery *> queriesWithPaths;
    for (int i = 0; i < queryCount; ++i)
    {
        BlastQuery * query = g_blastSearch->m_blastQueries.m_queries[i];
        if (query->getPathCount() == 0)
            out << "No query paths found for " << query->getName() << endl;
        else
            queriesWithPaths.push_back(query);
    }

    //Check that at least two queries have paths.
    int queriesWithPathsCount = queriesWithPaths.size();
    if (queriesWithPathsCount < 2)
    {
        err << "Bandage error: less than two queries have paths - no possible comparisons." << endl;
        return 1;
    }

    //Output a table header
    out << "Query 1 name\t";
    out << "Query 2 name\t";
    out << "Query 1 path\t";
    out << "Query 2 path\t";
    out << "Orientation\t";
    out << "Distance\t";
    out << "Distance path";
    if (sequences)
        out << "\tPath sequence";
    out << endl;

    //Look at each pair of queries.
    for (int i = 0; i < queriesWithPathsCount; ++i)
    {
        BlastQuery * query1 = queriesWithPaths[i];
        QString query1Name = query1->getName();

        QList<BlastQueryPath> query1Paths;
        if (allQueryPaths)
            query1Paths = query1->getPaths();
        else
            query1Paths = query1->getBestPaths();

        for (int j = i + 1; j < queriesWithPathsCount; ++j)
        {
            BlastQuery * query2 = queriesWithPaths[j];
            QString query2Name = query2->getName();

            QList<BlastQueryPath> query2Paths;
            if (allQueryPaths)
                query2Paths = query2->getPaths();
            else
                query2Paths = query2->getBestPaths();

            if (g_settings->findAllDistancePaths)
            {
                g_assemblyGraph->findAllPaths(query1Paths, query2Paths,
                                              g_settings->distanceOverlapOrientationSameStrand,
                                              g_settings->distanceOverlapOrientationOppositeStrand,
                                              g_settings->distanceOrientation1,
                                              g_settings->distanceOrientation2,
                                              g_settings->distanceOrientation3,
                                              g_settings->distanceOrientation4,
                                              g_settings->distancePathSearchDepth,
                                              g_settings->minDistancePathLength,
                                              g_settings->maxDistancePathLength);
            }
            else //Find shortest paths
            {
                g_assemblyGraph->findShortestPaths(query1Paths, query2Paths,
                                                   g_settings->distanceOverlapOrientationSameStrand,
                                                   g_settings->distanceOverlapOrientationOppositeStrand,
                                                   g_settings->distanceOrientation1,
                                                   g_settings->distanceOrientation2,
                                                   g_settings->distanceOrientation3,
                                                   g_settings->distanceOrientation4);
            }

            //Output the results
            int distancePathCount = g_memory->distanceSearchResults.size();
            for (int k = 0; k < distancePathCount; ++k)
            {
                out << query1Name << "\t";
                out << query2Name << "\t";
                out << g_memory->distanceSearchResults[k].m_query1Path.getString(true) << "\t";
                out << g_memory->distanceSearchResults[k].m_query2Path.getString(true) << "\t";
                out << g_memory->distanceSearchResults[k].m_orientation << "\t";
                out << g_memory->distanceSearchResults[k].m_distance << "\t";
                out << g_memory->distanceSearchResults[k].m_distancePath.getString(true);
                if (sequences)
                    out << "\t" << g_memory->distanceSearchResults[k].m_distancePath.getPathSequence();
                out << endl;
            }
        }
    }

    deleteBlastTempDirectory();
    return 0;
}



void printDistanceUsage(QTextStream * out, bool all)
{
    QStringList text;

    text << "Bandage distance takes two queries as input and will output (to stdout) the";
    text << "possible orientations and distances between them in the graph.";
    text << "";
    text << "Usage:    Bandage distance <graph> <queries> [options]";
    text << "";
    text << "Positional parameters:";
    text << "<graph>             A graph file of any type supported by Bandage";
    text << "<queries>           A FASTA file with at least two sequences (the distance search will be carried out for each pair of sequences in the file)";
    text << "";
    text << "Options:  --allquerypaths     Use all possible query paths in the graph for the distance search. If this option is not used, Bandage will use only the single best path for each query.";
    text << "--sequences         Include the path sequence in the output (note: can make output much larger)";
    text << "";
    text << "The options in the 'BLAST query paths' and 'Distance between queries' sections of the full Bandage settings (viewable via --helpall) are also relevant for the use of 'Bandage distance'.";
    text << "";

    getCommonHelp(&text);
    if (all)
        getSettingsUsage(&text);
    getOnlineHelpMessage(&text);

    outputText(text, out);
}



QString checkForInvalidDistanceOptions(QStringList arguments)
{
    checkOptionWithoutValue("--allquerypaths", &arguments);
    checkOptionWithoutValue("--sequences", &arguments);

    return checkForInvalidOrExcessSettings(&arguments);
}



void parseDistanceOptions(QStringList arguments, bool * allQueryPaths,
                          bool * sequences)
{
    *allQueryPaths = isOptionPresent("--allquerypaths", &arguments);
    *sequences = isOptionPresent("--sequences", &arguments);

    parseSettings(arguments);
}


