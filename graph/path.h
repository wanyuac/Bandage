//Copyright 2016 Ryan Wick

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


#ifndef PATH_H
#define PATH_H

#include <QByteArray>
#include <QList>
#include <vector>
#include <QString>
#include <QStringList>
#include "../program/globals.h"
#include "graphlocation.h"

class DeBruijnNode;
class DeBruijnEdge;

class Path
{
public:
    //CREATORS
    Path() {}
    static Path makeFromUnorderedNodes(QList<DeBruijnNode *> nodes,
                                       bool strandSpecific);
    static Path makeFromUnorderedNodes(std::vector<DeBruijnNode *> nodes,
                                       bool strandSpecific);
    static Path makeFromOrderedNodes(QList<DeBruijnNode *> nodes,
                                     bool circular);
    static Path makeFromString(QString pathString, bool circular,
                               QString * pathStringFailure);
    static Path makePathToEndOfNode(GraphLocation startLocation);

    //ACCESSORS
    QList<DeBruijnNode *> getNodes() const {return m_nodes;}
    QList<DeBruijnEdge *> getEdges() const {return m_edges;}
    bool isEmpty() const {return m_nodes.empty();}
    bool isCircular() const;
    bool haveSameNodes(Path other) const;
    bool hasNodeSubset(Path other) const;
    QByteArray getPathSequence() const;
    QString getFasta() const;
    QString getString(bool spaces) const;
    int getLength() const;
    QList<Path> extendPathInAllPossibleWays() const;
    bool canNodeFitOnEnd(DeBruijnNode * node, Path * extendedPath) const;
    bool canNodeFitAtStart(DeBruijnNode * node, Path * extendedPath) const;
    double getMeanDepth() const;
    bool containsNode(DeBruijnNode * node) const;
    bool containsEntireNode(DeBruijnNode * node) const;
    bool isInMiddleOfPath(DeBruijnNode * node) const;
    int numberOfOccurrencesInMiddleOfPath(DeBruijnNode * node) const;
    bool isStartingNode(DeBruijnNode * node) const;
    bool isEndingNode(DeBruijnNode * node) const;
    double getStartFraction() const;
    double getEndFraction() const;
    int getNodeCount() const;
    GraphLocation getStartLocation() const {return m_startLocation;}
    GraphLocation getEndLocation() const {return m_endLocation;}
    bool overlapsSameStrand(Path * other) const;
    bool overlapsOppositeStrand(Path * other) const;
    bool operator==(Path const &other) const;
    bool containsLocation(GraphLocation location) const;
    bool addNode(DeBruijnNode * newNode, bool strandSpecific, bool makeCircularIfPossible);
    void addNodeToEndNoCheck(DeBruijnNode * newNode);
    void extendPathToIncludeEntirityOfNodes();
    void setEndLocation(GraphLocation newEnd) {m_endLocation = newEnd;}
    void trimOneBaseOffEachEnd();

    //STATIC
    static QList<Path> getAllPossiblePaths(GraphLocation startLocation,
                                           GraphLocation endLocation,
                                           int nodeSearchDepth,
                                           int minDistance, int maxDistance,
                                           bool trimEnds = false);
    static QList<Path> getShortestPath(GraphLocation startLocation,
                                       GraphLocation endLocation);

private:
    GraphLocation m_startLocation;
    GraphLocation m_endLocation;
    QList<DeBruijnNode *> m_nodes;
    QList<DeBruijnEdge *> m_edges;

    void buildUnambiguousPathFromNodes(QList<DeBruijnNode *> nodes,
                                       bool strandSpecific);
    QByteArray modifySequenceUsingOverlap(QByteArray sequence, int overlap) const;
    bool checkForOtherEdges();
    void trimOneBaseOffStart();
    void trimOneBaseOffEnd();
};

#endif // PATH_H
