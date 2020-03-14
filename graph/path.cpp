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


#include "path.h"
#include "debruijnnode.h"
#include "debruijnedge.h"
#include "../blast/blasthit.h"
#include "../blast/blastquery.h"
#include <QRegularExpression>
#include "assemblygraph.h"
#include <QStringList>
#include <QApplication>
#include <QHash>
#include <QSet>
#include <set>
#include <utility>
#include <limits>


//These will try to produce a path using the given nodes.
//They will only succeed if the nodes produce one and only one path.
//If they are disconnected, branching or ambiguous, they will fail
//and only contruct an empty Path.
Path Path::makeFromUnorderedNodes(QList<DeBruijnNode *> nodes,
                                  bool strandSpecific)
{
    Path path;
    path.buildUnambiguousPathFromNodes(nodes, strandSpecific);
    return path;
}

Path Path::makeFromUnorderedNodes(std::vector<DeBruijnNode *> nodes,
                                  bool strandSpecific)
{
    QList<DeBruijnNode *> nodesList;
    for (size_t i = 0; i < nodes.size(); ++i)
        nodesList.push_back(nodes[i]);

    Path path;
    path.buildUnambiguousPathFromNodes(nodesList, strandSpecific);
    return path;
}


//This will build a Path from an ordered list of nodes.  If the nodes
//form a valid path (i.e. there is an edge connecting each step along
//the way), a Path is made, otherwise just an empty Path is made.
//This function needs exact, strand-specific nodes.  If circular is
//given, then it will also look for an edge connecting the last node
//to the first.
Path Path::makeFromOrderedNodes(QList<DeBruijnNode *> nodes, bool circular)
{
    Path path;

    path.m_nodes = nodes;

    int targetNumberOfEdges = path.m_nodes.size() - 1;
    if (circular)
        ++targetNumberOfEdges;

    for (int i = 0; i < targetNumberOfEdges; ++i)
    {
        int firstNodeIndex = i;
        int secondNodeIndex = i + 1;
        if (secondNodeIndex >= path.m_nodes.size())
            secondNodeIndex = 0;

        DeBruijnNode * node1 = path.m_nodes[firstNodeIndex];
        DeBruijnNode * node2 = path.m_nodes[secondNodeIndex];

        bool foundEdge = false;
        const std::vector<DeBruijnEdge *> * edges = node1->getEdgesPointer();
        for (size_t j = 0; j < edges->size(); ++j)
        {
            DeBruijnEdge * edge = (*edges)[j];
            if (edge->getStartingNode() == node1 && edge->getEndingNode() == node2)
            {
                path.m_edges.push_back(edge);
                foundEdge = true;
                break;
            }
        }

        //If we failed to find an edge connecting the nodes, then
        //the path failed.
        if (!foundEdge)
        {
            path.m_nodes.clear();
            path.m_edges.clear();
            return path;
        }
    }

    if (path.m_nodes.empty())
        return path;

    //If the code got here, then the path building was successful.
    path.m_startLocation = GraphLocation::startOfNode(path.m_nodes.front());
    path.m_endLocation = GraphLocation::endOfNode(path.m_nodes.back());

    return path;
}



Path Path::makeFromString(QString pathString, bool circular,
                          QString * pathStringFailure)
{
    Path path;

    QRegularExpression re("^(?:\\(([0-9]+)\\) ?)*((?:[^,]+[-\\+], ?)*[^,]+[-\\+])(?: ?\\(([0-9]+)\\))*$");
    QRegularExpressionMatch match = re.match(pathString);

    //If the string failed to match the regex, return an empty path.
    if (!match.hasMatch())
    {
        *pathStringFailure = "the text is not formatted correctly";
        return path;
    }

    QString startPosString = match.captured(1);
    QString nodeListString = match.captured(2);
    QString endPosString = match.captured(3);

    //Circular paths cannot have start and end positions.
    if (circular && (startPosString != "" || endPosString != ""))
    {
        *pathStringFailure = "circular paths cannot contain start or end positions";
        return path;
    }

    //Make sure there is at least one proposed node name listed.
    QStringList nodeNameList = nodeListString.simplified().split(",", QString::SkipEmptyParts);
    if (nodeNameList.empty())
    {
        *pathStringFailure = "the text is not formatted correctly";
        return path;
    }

    //Find which node names are and are not actually in the graph. 
    QList<DeBruijnNode *> nodesInGraph;
    QStringList nodesNotInGraph;
    for (int i = 0; i < nodeNameList.size(); ++i)
    {
        QString nodeName = nodeNameList[i].simplified();
        if (g_assemblyGraph->m_deBruijnGraphNodes.contains(nodeName))
            nodesInGraph.push_back(g_assemblyGraph->m_deBruijnGraphNodes[nodeName]);
        else
            nodesNotInGraph.push_back(nodeName);
    }

    //If the path contains nodes not in the graph, we fail.
    if (nodesNotInGraph.size() > 0)
    {
        *pathStringFailure = "the following nodes are not in the graph: ";
        for (int i = 0; i < nodesNotInGraph.size(); ++i)
        {
            *pathStringFailure += nodesNotInGraph[i];
            if (i != nodesNotInGraph.size() - 1)
                *pathStringFailure += ", ";
        }
        return path;
    }

    //If the code got here, then the list at least consists of valid nodes.
    //We now use it to create a Path object.
    path = Path::makeFromOrderedNodes(nodesInGraph, circular);

    //If the path is empty, then we don't have to worry about start/end
    //positions and we just return it.
    if (path.isEmpty())
    {
        if (circular)
            *pathStringFailure = "the nodes do not form a circular path";
        else
            *pathStringFailure = "the nodes do not form a path";
        return path;
    }

    //If the code got here, then a path was made, and now we must check whether
    //the start/end points are valid.
    DeBruijnNode * firstNode = path.m_nodes.front();
    DeBruijnNode * lastNode = path.m_nodes.back();

    if (startPosString.length() > 0)
    {
        int startPos = startPosString.toInt();
        if (startPos < 1 || startPos > firstNode->getLength())
        {
            *pathStringFailure = "starting node position not valid";
            return Path();
        }

        path.m_startLocation = GraphLocation(firstNode, startPos);
    }
    else
        path.m_startLocation = GraphLocation::startOfNode(firstNode);


    if (endPosString.length() > 0)
    {
        int endPos = endPosString.toInt();
        if (endPos < 1 || endPos > lastNode->getLength())
        {
            *pathStringFailure = "ending node position not valid";
            return Path();
        }

        path.m_endLocation = GraphLocation(lastNode, endPos);
    }
    else
        path.m_endLocation = GraphLocation::endOfNode(lastNode);

    return path;
}



Path Path::makePathToEndOfNode(GraphLocation startLocation)
{
    Path path;

    DeBruijnNode * node = startLocation.getNode();
    path.m_nodes.push_back(node);
    path.m_startLocation = startLocation;
    path.m_endLocation = GraphLocation::endOfNode(node);

    return path;
}



void Path::buildUnambiguousPathFromNodes(QList<DeBruijnNode *> nodes,
                                         bool strandSpecific)
{
    if (nodes.isEmpty())
        return;

    //Loop through the nodes, trying to add them to the Path.  If a node can't
    //be added, then we fail and make an empty Path.  If one can be added, we
    //quit the loop and try again with the remaining nodes.
    while (nodes.size() > 0)
    {
        bool addSuccess = false;
        for (int i = 0; i < nodes.size(); ++i)
        {
            addSuccess = addNode(nodes.at(i), strandSpecific, true);
            if (addSuccess)
            {
                nodes.removeAt(i);
                break;
            }
        }

        if (!addSuccess)
        {
            m_nodes.clear();
            m_edges.clear();
            return;
        }
    }

    //If the nodes in the path contain other edges which connect them to each
    //other, then the path is ambiguous and we fail.
    if (checkForOtherEdges())
    {
        m_nodes.clear();
        m_edges.clear();
        return;
    }

    if (m_nodes.empty())
        return;

    //If the code got here, then the path building was successful.
    m_startLocation = GraphLocation::startOfNode(m_nodes.front());
    m_endLocation = GraphLocation::endOfNode(m_nodes.back());
}


//This function will try to add a node to the path on either end.
//It will only succeed (and return true) if there is a single way
//to add the node on one of the path's ends.
//It can, however, add a node that connects the end to both ends,
//making a circular Path.
bool Path::addNode(DeBruijnNode * newNode, bool strandSpecific, bool makeCircularIfPossible)
{
    //If the Path is empty, then this function always succeeds.
    if (m_nodes.isEmpty())
    {
        m_nodes.push_back(newNode);
        m_startLocation = GraphLocation::startOfNode(newNode);
        m_endLocation = GraphLocation::endOfNode(newNode);

        if (makeCircularIfPossible)
        {
            //If there is an edge connecting the node to itself, then add that
            //too to make a circular path.
            DeBruijnEdge * selfLoopingEdge = newNode->getSelfLoopingEdge();
            if (selfLoopingEdge != 0)
                m_edges.push_back(selfLoopingEdge);
        }

        return true;
    }

    //If the Path is circular, then this function fails, as there
    //is no way to add a node to a circular path without making
    //it ambiguous.
    if (isCircular())
        return false;

    //Check to see if the node can be added anywhere in the middle
    //of the Path.  If so, this function fails.
    for (int i = 1; i < m_nodes.size() - 1; ++i)
    {
        DeBruijnNode * middleNode = m_nodes.at(i);
        if (middleNode->isNodeConnected(newNode))
            return false;
    }

    DeBruijnNode * firstNode = m_nodes.front();
    DeBruijnNode * lastNode = m_nodes.back();

    DeBruijnEdge * edgeIntoFirst = firstNode->doesNodeLeadIn(newNode);
    DeBruijnEdge * edgeAwayFromLast = lastNode->doesNodeLeadAway(newNode);

    //If not strand-specific, then we also check to see if the reverse
    //complement of the new node can be added.
    DeBruijnEdge * revCompEdgeIntoFirst = 0;
    DeBruijnEdge * revCompEdgeAwayFromLast = 0;
    if (!strandSpecific)
    {
        revCompEdgeIntoFirst = firstNode->doesNodeLeadIn(newNode->getReverseComplement());
        revCompEdgeAwayFromLast = lastNode->doesNodeLeadAway(newNode->getReverseComplement());
    }

    //To be successful, either:
    // 1) exactly one of the four edge pointers should be non-null.  This
    //    indicates the node extends a linear path.
    // 2) there is both an edge away from the last and an edge into the first.
    //    This indicates that the node completes a circular Path.
    if (edgeIntoFirst == 0 && edgeAwayFromLast == 0 &&
            revCompEdgeIntoFirst == 0 && revCompEdgeAwayFromLast == 0)
        return false;

    if (edgeIntoFirst != 0 && edgeAwayFromLast == 0 &&
            revCompEdgeIntoFirst == 0 && revCompEdgeAwayFromLast == 0)
    {
        m_nodes.push_front(newNode);
        m_startLocation = GraphLocation::startOfNode(newNode);
        m_edges.push_front(edgeIntoFirst);
        return true;
    }

    if (edgeIntoFirst == 0 && edgeAwayFromLast != 0 &&
            revCompEdgeIntoFirst == 0 && revCompEdgeAwayFromLast == 0)
    {
        m_nodes.push_back(newNode);
        m_endLocation = GraphLocation::endOfNode(newNode);
        m_edges.push_back(edgeAwayFromLast);
        return true;
    }

    if (edgeIntoFirst == 0 && edgeAwayFromLast == 0 &&
            revCompEdgeIntoFirst != 0 && revCompEdgeAwayFromLast == 0)
    {
        newNode = newNode->getReverseComplement();
        m_nodes.push_front(newNode);
        m_startLocation = GraphLocation::startOfNode(newNode);
        m_edges.push_front(revCompEdgeIntoFirst);
        return true;
    }

    if (edgeIntoFirst == 0 && edgeAwayFromLast == 0 &&
            revCompEdgeIntoFirst == 0 && revCompEdgeAwayFromLast != 0)
    {
        newNode = newNode->getReverseComplement();
        m_nodes.push_back(newNode);
        m_endLocation = GraphLocation::endOfNode(newNode);
        m_edges.push_back(revCompEdgeAwayFromLast);
        return true;
    }

    if (edgeIntoFirst != 0 && edgeAwayFromLast != 0 &&
            revCompEdgeIntoFirst == 0 && revCompEdgeAwayFromLast == 0)
    {
        m_edges.push_back(edgeAwayFromLast);
        m_nodes.push_back(newNode);
        m_edges.push_back(edgeIntoFirst);
        return true;
    }

    if (edgeIntoFirst == 0 && edgeAwayFromLast == 0 &&
            revCompEdgeIntoFirst != 0 && revCompEdgeAwayFromLast != 0)
    {
        m_edges.push_back(revCompEdgeAwayFromLast);
        m_nodes.push_back(newNode->getReverseComplement());
        m_edges.push_back(revCompEdgeIntoFirst);
        return true;
    }

    //If the code got here, then there are multiple ways of adding the node, so
    //we fail.
    return false;
}



//This function adds a node to the end of the path.  It assumes that the node
//can in fact be added to the end - it does not check this.
void Path::addNodeToEndNoCheck(DeBruijnNode * newNode)
{
    m_nodes.push_back(newNode);
    if (m_nodes.size() > 1)
    {
        DeBruijnNode * secondToLastNode = m_nodes[m_nodes.size() - 2];
        m_edges.push_back(secondToLastNode->doesNodeLeadAway(newNode));
    }
    m_endLocation = GraphLocation::endOfNode(newNode);
}


//This function looks to see if there are other edges connecting path nodes
//that aren't in the list of path edges.  If so, it returns true.
//This is used to check whether a Path is ambiguous or node.
bool Path::checkForOtherEdges()
{
    //First compile a list of all edges which connect any
    //node in the Path to any other node in the Path.
    QList<DeBruijnEdge *> allConnectingEdges;
    for (int i = 0; i < m_nodes.size(); ++i)
    {
        DeBruijnNode * startingNode = m_nodes[i];
        const std::vector<DeBruijnEdge *> * startingNodeEdges = startingNode->getEdgesPointer();
        for (int j = 0; j < m_nodes.size(); ++j)
        {
            DeBruijnNode * endingNode = m_nodes[j];
            for (size_t k = 0; k < startingNodeEdges->size(); ++k)
            {
                DeBruijnEdge * edge = (*startingNodeEdges)[k];
                if (edge->getStartingNode() == startingNode &&
                        edge->getEndingNode() == endingNode)
                    allConnectingEdges.push_back(edge);
            }
        }
    }

    //If this list of edges is greater than the number of edges in the path,
    //then other edges exist.
    return allConnectingEdges.size() > m_edges.size();
}




//This function extracts the sequence for the whole path.  It uses the overlap
//value in the edges to remove sequences that are duplicated at the end of one
//node and the start of the next.
QByteArray Path::getPathSequence() const
{
    if (m_nodes.empty())
        return "";

    QByteArray sequence;
    QByteArray firstNodeSequence = m_nodes[0]->getSequence();

    //If the path is circular, we trim the overlap from the first node.
    if (isCircular())
    {
        int overlap = m_edges.back()->getOverlap();
        if (overlap != 0)
            firstNodeSequence = modifySequenceUsingOverlap(firstNodeSequence, overlap);
        sequence += firstNodeSequence;
    }

    //If the path is linear, then we begin either with the entire first node
    //sequence or part of it.
    else
    {
        int rightChars = firstNodeSequence.length() - m_startLocation.getPosition() + 1;
        sequence += firstNodeSequence.right(rightChars);
    }

    //The middle nodes are not affected by whether or not the path is circular
    //or has partial node ends.
    for (int i = 1; i < m_nodes.size(); ++i)
    {
        int overlap = m_edges[i-1]->getOverlap();
        QByteArray nodeSequence = m_nodes[i]->getSequence();
        if (overlap != 0)
            nodeSequence = modifySequenceUsingOverlap(nodeSequence, overlap);
        sequence += nodeSequence;
    }

    DeBruijnNode * lastNode = m_nodes.back();
    int amountToTrimFromEnd = lastNode->getLength() - m_endLocation.getPosition();
    sequence.chop(amountToTrimFromEnd);

    return sequence;
}


//This function will trim bases from the start of a sequence (in the case of
//positive overlap) or add Ns to the start (in the case of negative overlap).
QByteArray Path::modifySequenceUsingOverlap(QByteArray sequence, int overlap) const
{
    if (overlap > 0)
    {
        int rightChars = sequence.length() - overlap;
        if (rightChars >= 0)
            sequence = sequence.right(rightChars);
    }
    else if (overlap < 0)
        sequence = QByteArray(-overlap, 'N') + sequence;

    return sequence;
}


int Path::getLength() const
{
    if (m_nodes.empty())
        return 0;

    int length = 0;
    for (int i = 0; i < m_nodes.size(); ++i)
        length += m_nodes[i]->getLength();

    for (int i = 0; i < m_edges.size(); ++i)
        length -= m_edges[i]->getOverlap();

    length -= m_startLocation.getPosition() - 1;

    DeBruijnNode * lastNode = m_nodes.back();
    length -= lastNode->getLength() - m_endLocation.getPosition();

    return length;
}


QString Path::getFasta() const
{
    //The description line is a comma-delimited list of the nodes in the path
    QString fasta = ">" + getString(false);

    if (isCircular())
        fasta += "(circular)";
    fasta += "\n";

    QString pathSequence = getPathSequence();
    int charactersOnLine = 0;
    for (int i = 0; i < pathSequence.length(); ++i)
    {
        fasta += pathSequence.at(i);
        ++charactersOnLine;
        if (charactersOnLine >= 70)
        {
            fasta += "\n";
            charactersOnLine = 0;
        }
    }
    fasta += "\n";

    return fasta;
}



QString Path::getString(bool spaces) const
{
    QString output;
    for (int i = 0; i < m_nodes.size(); ++i)
    {
        if (i == 0 && !m_startLocation.isAtStartOfNode())
        {
            output += "(" + QString::number(m_startLocation.getPosition()) + ")";
            if (spaces)
                output += " ";
        }

        output += m_nodes[i]->getName();
        if (i < m_nodes.size() - 1)
        {
            output += ",";
            if (spaces)
                output += " ";
        }

        if (i == m_nodes.size() - 1 && !m_endLocation.isAtEndOfNode())
        {
            if (spaces)
                output += " ";
            output += "(" + QString::number(m_endLocation.getPosition()) + ")";
        }
    }
    return output;
}


//This function tests whether the last node in the path leads into the first.
bool Path::isCircular() const
{
    if (isEmpty())
        return false;
    if (m_edges.empty())
        return false;

    //A circular path should have the same number of nodes and edges.
    if (m_nodes.size() != m_edges.size())
        return false;

    DeBruijnEdge * lastEdge = m_edges.back();
    DeBruijnNode * firstNode = m_nodes.front();
    DeBruijnNode * lastNode = m_nodes.back();

    return (lastEdge->getStartingNode() == lastNode &&
            lastEdge->getEndingNode() == firstNode);
}



//These functions test whether the specified node could be added to
//the end/front of the path to form a larger valid path.
//If so, they set the path pointed to by extendedPath to equal the new, larger
//path.
bool Path::canNodeFitOnEnd(DeBruijnNode * node, Path * extendedPath) const
{
    if (isEmpty())
    {
        QList<DeBruijnNode *> nodeList;
        nodeList.push_back(node);
        *extendedPath = Path::makeFromOrderedNodes(nodeList, false);
        return true;
    }
    if (isCircular())
        return false;

    DeBruijnNode * lastNode = m_nodes.back();
    const std::vector<DeBruijnEdge *> * lastNodeEdges = lastNode->getEdgesPointer();
    for (size_t i = 0; i < lastNodeEdges->size(); ++i)
    {
        DeBruijnEdge * edge = (*lastNodeEdges)[i];
        if (edge->getStartingNode() == lastNode && edge->getEndingNode() == node)
        {
            *extendedPath = *this;
            extendedPath->m_edges.push_back(edge);
            extendedPath->m_nodes.push_back(node);
            extendedPath->m_endLocation = GraphLocation::endOfNode(node);
            return true;
        }
    }

    return false;
}

bool Path::canNodeFitAtStart(DeBruijnNode * node, Path * extendedPath) const
{
    if (isEmpty())
    {
        QList<DeBruijnNode *> nodeList;
        nodeList.push_back(node);
        *extendedPath = Path::makeFromOrderedNodes(nodeList, false);
        return true;
    }
    if (isCircular())
        return false;

    DeBruijnNode * firstNode = m_nodes.front();
    const std::vector<DeBruijnEdge *> * firstNodeEdges = firstNode->getEdgesPointer();
    for (size_t i = 0; i < firstNodeEdges->size(); ++i)
    {
        DeBruijnEdge * edge = (*firstNodeEdges)[i];
        if (edge->getStartingNode() == node && edge->getEndingNode() == firstNode)
        {
            *extendedPath = *this;
            extendedPath->m_edges.push_front(edge);
            extendedPath->m_nodes.push_front(node);
            extendedPath->m_startLocation = GraphLocation::startOfNode(node);
            return true;
        }
    }

    return false;
}


//This function builds all possible paths between the given start and end,
//within the given restrictions.
QList<Path> Path::getAllPossiblePaths(GraphLocation startLocation,
                                      GraphLocation endLocation,
                                      int nodeSearchDepth,
                                      int minDistance, int maxDistance,
                                      bool trimEnds)
{
    QList<Path> finishedPaths;
    QList<Path> unfinishedPaths;

    Path startPath;
    startPath.addNode(startLocation.getNode(), true, false);
    startPath.m_startLocation = startLocation;
    startPath.m_endLocation = GraphLocation::endOfNode(startLocation.getNode());
    unfinishedPaths.push_back(startPath);

    for (int i = 0; i <= nodeSearchDepth; ++i)
    {
        QApplication::processEvents();

        //Look at each of the unfinished paths to see if they end with the end
        //node.  If so, see if it has the appropriate length.
        //If it does, it will go into the final returned list.
        //If it doesn't and it's over length, then it will be removed.
        QList<Path>::iterator j = unfinishedPaths.begin();
        while (j != unfinishedPaths.end())
        {
            DeBruijnNode * lastNode = (*j).m_nodes.back();
            if (lastNode == endLocation.getNode() &&
                    ((*j).m_nodes.size() > 1 || endLocation.getPosition() > startLocation.getPosition()))
            {
                Path potentialFinishedPath = *j;
                potentialFinishedPath.m_endLocation = endLocation;
                if (trimEnds)
                    potentialFinishedPath.trimOneBaseOffEachEnd();
                int length = potentialFinishedPath.getLength();
                if (length >= minDistance && length <= maxDistance)
                    finishedPaths.push_back(potentialFinishedPath);
                ++j;
            }
            else
            {
                if ((*j).getLength() > maxDistance)
                    j = unfinishedPaths.erase(j);
                else
                    ++j;
            }
        }

        //Make new unfinished paths by extending each of the paths.
        QList<Path> newUnfinishedPaths;
        for (int j = 0; j < unfinishedPaths.size(); ++j)
            newUnfinishedPaths.append(unfinishedPaths[j].extendPathInAllPossibleWays());
        unfinishedPaths = newUnfinishedPaths;
    }

    return finishedPaths;
}


//This function takes the current path and extends it in all possible ways by
//adding one more node, then returning a list of the new paths.  How many paths
//it returns depends on the number of edges leaving the last node in the path.
QList<Path> Path::extendPathInAllPossibleWays() const
{
    QList<Path> returnList;

    if (isEmpty())
        return returnList;

    //Since circular paths don't really have an end to extend, this function
    //doesn't work for them.
    if (isCircular())
        return returnList;

    DeBruijnNode * lastNode = m_nodes.back();
    std::vector<DeBruijnEdge *> nextEdges = lastNode->getLeavingEdges();
    for (size_t i = 0; i < nextEdges.size(); ++i)
    {
        DeBruijnEdge * nextEdge = nextEdges[i];
        DeBruijnNode * nextNode = nextEdge->getEndingNode();

        Path newPath(*this);
        newPath.m_edges.push_back(nextEdge);
        newPath.m_nodes.push_back(nextNode);
        newPath.m_endLocation = GraphLocation::endOfNode(nextNode);

        returnList.push_back(newPath);
    }

    return returnList;
}



double Path::getMeanDepth() const
{
    long double depthTimesLengthSum = 0.0;
    int nodeLengthTotal = 0;
    for (int i = 0; i < m_nodes.size(); ++i)
    {
        DeBruijnNode * node = m_nodes[i];
        depthTimesLengthSum += node->getDepth() * node->getLength();
        nodeLengthTotal += node->getLength();
    }

    return depthTimesLengthSum / nodeLengthTotal;
}



bool Path::operator==(Path const &other) const
{
    return (m_nodes == other.m_nodes &&
            m_startLocation == other.m_startLocation &&
            m_endLocation == other.m_endLocation);
}


bool Path::haveSameNodes(Path other) const
{
    return (m_nodes == other.m_nodes);
}



//This function checks to see if this path is actually a sub-path (i.e.
//entirely contained within) the other given path.
//It ignores start/end type and position, looking only at the nodes.
//If the two paths have the same nodes, it will return false.
bool Path::hasNodeSubset(Path other) const
{
    //To contain this path, the other path should be have a larger number of
    //nodes.
    int nodeCountDifference = other.m_nodes.size() - m_nodes.size();
    if (nodeCountDifference <= 0)
        return false;
    
    //If the paths have the same number of nodes, check to see if they are
    //identical.
    if (nodeCountDifference == 0)
        return (m_nodes == other.m_nodes);
    
    //If the code got here, then the other path has more nodes than this path.
    //We now see if we can find an ordered set of nodes in the other path that
    //matches this path's nodes.
    for (int i = 0; i <= nodeCountDifference; ++i)
    {
        QList<DeBruijnNode *> otherPathNodeSubset = other.m_nodes.mid(i, m_nodes.size());
        if (m_nodes == otherPathNodeSubset)
            return true;
    }
    
    return false;
}



void Path::extendPathToIncludeEntirityOfNodes()
{
    if (m_nodes.empty())
        return;

    m_startLocation = GraphLocation::startOfNode(m_nodes.front());
    m_endLocation = GraphLocation::endOfNode(m_nodes.back());
}


bool Path::containsNode(DeBruijnNode * node) const
{
    return m_nodes.contains(node);
}

bool Path::containsEntireNode(DeBruijnNode * node) const
{
    if (m_nodes.empty())
        return false;

    if (m_nodes.size() == 1) {
        if (m_nodes.front() != node)
            return false;
        return m_startLocation.isAtStartOfNode() && m_endLocation.isAtEndOfNode();
    }

    if (m_nodes.front() == node && m_startLocation.isAtStartOfNode())
        return true;

    if (m_nodes.back() == node && m_endLocation.isAtEndOfNode())
        return true;

    for (int i = 1; i < m_nodes.size() - 1; ++i)
    {
        if (m_nodes[i] == node)
            return true;
    }

    return false;
}




bool Path::isInMiddleOfPath(DeBruijnNode * node) const
{
    return containsNode(node) && !isStartingNode(node) && !isEndingNode(node);
}


//This function counts the number of times the node is in the path, not
//counting the first or last nodes.
int Path::numberOfOccurrencesInMiddleOfPath(DeBruijnNode * node) const
{
    int count = 0;
    for (int i = 1; i < m_nodes.size() - 1; ++i)
    {
        if (m_nodes[i] == node)
            ++count;
    }
    return count;
}

bool Path::isStartingNode(DeBruijnNode * node) const
{
    if (m_nodes.empty())
        return false;
    return m_nodes.front() == node;
}

bool Path::isEndingNode(DeBruijnNode * node) const
{
    if (m_nodes.empty())
        return false;
    return m_nodes.back() == node;
}


double Path::getStartFraction() const
{
    if (m_nodes.empty())
        return 0.0;

    int firstNodeLength = m_nodes.front()->getLength();
    if (firstNodeLength == 0)
        return 0.0;

    return double(m_startLocation.getPosition() - 1) / firstNodeLength;
}

double Path::getEndFraction() const
{
    if (m_nodes.empty())
        return 1.0;

    int lastNodeLength = m_nodes.back()->getLength();
    if (lastNodeLength == 0)
        return 1.0;

    return double(m_endLocation.getPosition()) / lastNodeLength;
}


int Path::getNodeCount() const
{
    return m_nodes.size();
}



//This function uses an implementation of Dijkstra's algorithm to find the
//shortest path between two locations.
//It returns a list of length zero (if it failed to find a shortest path) or
//a list of length one (if it succeeded).
QList<Path> Path::getShortestPath(GraphLocation startLocation,
                                  GraphLocation endLocation)
{
    DeBruijnNode * startingNode = startLocation.getNode();
    DeBruijnNode * destinationNode = endLocation.getNode();

    const int INF = std::numeric_limits<int>::max();

    //We need to handle the special case of when our destination is before our
    //start but in the same node.  In this situation, we must visit the first
    //node twice.
    bool visitFirstNodeTwice = false;
    if (startingNode == destinationNode &&
            startLocation.getPosition() >= endLocation.getPosition())
        visitFirstNodeTwice = true;

    //Initialise the distances to infinity, initialise the paths to an empty
    //path and label all nodes as unvisited.
    QHash<DeBruijnNode *, int> distances;
    QHash<DeBruijnNode *, Path> paths;
    QMapIterator<QString, DeBruijnNode*> i(g_assemblyGraph->m_deBruijnGraphNodes);
    while (i.hasNext())
    {
        i.next();
        DeBruijnNode * node = i.value();

        distances.insert(node, INF);
        paths.insert(node, Path());
    }

    //Set the starting node distance to a non-INF value and set its path.
    int distanceFromStartLocationToEndOfStartNode = startingNode->getLength() - startLocation.getPosition();
    distances.insert(startingNode, distanceFromStartLocationToEndOfStartNode);
    paths.insert(startingNode, Path::makePathToEndOfNode(startLocation));

    //Create a priority queue of the nodes using their distance.
    std::set<std::pair<int, DeBruijnNode *> > nodeQueue;
    nodeQueue.insert(std::make_pair(distances.value(startingNode), startingNode));

    //The main Dijkstra's algorithm loop.
    while (!nodeQueue.empty())
    {
        //Note the current distance, node and path. The set is ordered, so this
        //will always be the node with the smallest distance.
        int distance = nodeQueue.begin()->first;
        DeBruijnNode * node = nodeQueue.begin()->second;
        Path path = paths.value(node);

        //If the destination node has been visited, then the search is
        //finished!
        if (node == destinationNode && !visitFirstNodeTwice)
        {
            Path shortestPath = path;
            shortestPath.setEndLocation(endLocation);
            shortestPath.trimOneBaseOffEachEnd();
            QList<Path> returnList;
            returnList.push_back(shortestPath);
            return returnList;
        }

        //Remove the current node from the queue.
        nodeQueue.erase(nodeQueue.begin());

        //If our start node is our end node (but with a later start position
        //than end position) and this is our first time visiting it, then we
        //need to reset its distance so it will get visited again.
        if (node == destinationNode && visitFirstNodeTwice)
        {
            distances.insert(node, INF);
            visitFirstNodeTwice = false;
        }

        //Look at each of the current node's neighbours.
        std::vector<DeBruijnNode *> neighbours = node->getDownstreamNodes();
        for (size_t i = 0; i < neighbours.size(); ++i)
        {
            DeBruijnNode * neighbour = neighbours[i];
            int oldDistance = distances.value(neighbour);
            int newDistance = distance + neighbour->getLength();

            //If the distance is an improvement...
            if (newDistance < oldDistance)
            {
                //Update the path for the neighbour
                Path newPath = path;
                newPath.addNodeToEndNoCheck(neighbour);
                paths.insert(neighbour, newPath);

                //Update the distance for the neighbour
                distances.insert(neighbour, newDistance);

                //The neighbour may or may not already be in the queue, but we
                //want to erase it if it is.
                nodeQueue.erase(std::make_pair(oldDistance, neighbour));

                //Add the neighbour to the queue, with the new distance.
                nodeQueue.insert(std::make_pair(newDistance, neighbour));
            }
        }
    }

    //If the code got here, then no shortest path was found, meaning that the
    //two locations are disconnected.  Return an empty path.
    return QList<Path>();
}



//This function checks to see if any part of the other path overlaps with this
//one and returns true if so.  It only checks the nodes in the paths, so it is
//a same-strand-only overlap check.
bool Path::overlapsSameStrand(Path * other) const
{
    QList<DeBruijnNode *> path1Nodes = getNodes();
    QList<DeBruijnNode *> path2Nodes = other->getNodes();

    //If either path is empty, then they can't overlap.
    if (path1Nodes.empty() || path2Nodes.empty())
        return false;

    bool path1Circular = isCircular();
    bool path2Circular = other->isCircular();


    //If any of the whole nodes in the first path are also whole nodes in the
    //second path, then the paths overlap.
    for (int i = 0; i < path1Nodes.size(); ++i)
    {
        DeBruijnNode * path1Node = path1Nodes[i];
        int path1NodeStart = 1;
        int path1NodeEnd = path1Node->getLength();

        //If this is the first and/or last node, then we update the start and/or
        //end positions.
        if (!path1Circular)
        {
            if (i == 0)
                path1NodeStart = getStartLocation().getPosition();
            if (i == path1Nodes.size() - 1)
                path1NodeEnd = getEndLocation().getPosition();
        }

        for (int j = 0; j < path2Nodes.size(); ++j)
        {
            DeBruijnNode * path2Node = path2Nodes[j];
            int path2NodeStart = 1;
            int path2NodeEnd = path2Node->getLength();

            //If this is the first and/or last node, then we update the start and/or
            //end positions.
            if (!path2Circular)
            {
                if (j == 0)
                    path2NodeStart = other->getStartLocation().getPosition();
                if (j == path2Nodes.size() - 1)
                    path2NodeEnd = other->getEndLocation().getPosition();
            }

            //If we find that the a node is in both paths, then we must check
            //to see if they overlap.
            if (path1Node == path2Node)
            {
                if (path1NodeEnd >= path2NodeStart && path2NodeEnd >= path1NodeStart)
                    return true;
            }
        }
    }

    //If the code got here, that means that none of the nodes in the paths
    //overlap.
    return false;
}


//This node checks whether two paths overlap on opposite strands.  This is a
//more difficult test, because for LastGraph graphs, the nodes are shifted, so
//there isn't an easy way to change a path into its reverse complement.
//It does the somewhat laborious task of taking each position in the path,
//grabbing the reverse complement position and checking to see if that's in the
//other path.
//Note that this isn't bullet-proof, as there can be multiple possible reverse
//complement positions in a LastGraph graph, but I think it is good enough.
bool Path::overlapsOppositeStrand(Path * other) const
{
    //If there is only one node in the path, then things are simpler.  We just
    //check each of the path positions in that node.
    if (m_nodes.size() == 1)
    {
        DeBruijnNode * node = m_startLocation.getNode();
        int startPosition = m_startLocation.getPosition();
        int endPosition = m_endLocation.getPosition();

        for (int i = startPosition; i <= endPosition; ++i)
        {
            GraphLocation pathLocation(node, i);
            GraphLocation reverseComplement = pathLocation.reverseComplementLocation();
            if (other->containsLocation(reverseComplement))
                return true;
        }

        return false;
    }

    //If there are multiple nodes in the path, then we must treat the first,
    //middle and last node differently.

    DeBruijnNode * startNode = m_startLocation.getNode();
    int startPosition = m_startLocation.getPosition();
    for (int i = startPosition; i <= startNode->getLength(); ++i)
    {
        GraphLocation pathLocation(startNode, i);
        GraphLocation reverseComplement = pathLocation.reverseComplementLocation();
        if (other->containsLocation(reverseComplement))
            return true;
    }

    for (int i = 1; i < m_nodes.size() - 1; ++i)
    {
        DeBruijnNode * middleNode = m_nodes[i];
        for (int i = 1; i <= middleNode->getLength(); ++i)
        {
            GraphLocation pathLocation(middleNode, i);
            GraphLocation reverseComplement = pathLocation.reverseComplementLocation();
            if (other->containsLocation(reverseComplement))
                return true;
        }
    }

    DeBruijnNode * endNode = m_endLocation.getNode();
    int endPosition = m_endLocation.getPosition();
    for (int i = 1; i <= endPosition; ++i)
    {
        GraphLocation pathLocation(endNode, i);
        GraphLocation reverseComplement = pathLocation.reverseComplementLocation();
        if (other->containsLocation(reverseComplement))
            return true;
    }

    return false;
}



//This function reduces a path size by two, by trimming a base off each end.  If
//the path only has 1 or 2 bases, this will result in an empty path.
void Path::trimOneBaseOffEachEnd()
{
    //Can't trim bases from circular paths.
    if (isCircular())
        return;

    trimOneBaseOffStart();
    trimOneBaseOffEnd();

    //Check to see if the path still has a length of at least 1.  If not, make
    //it an empty path.
    if (getLength() < 1)
    {
        m_nodes.clear();
        m_edges.clear();
    }
}


void Path::trimOneBaseOffStart()
{
    DeBruijnNode * startNode = m_startLocation.getNode();
    int startPosition = m_startLocation.getPosition();

    //If there is room in the first node to move forward by one, do that.
    if (startPosition < startNode->getLength())
        m_startLocation.moveLocation(1);

    //If not, we have to remove the start node and set the Path start to the
    //beginning of the next node.
    else
    {
        if (m_nodes.size() > 0)
            m_nodes.pop_front();
        if (m_edges.size() > 0)
            m_edges.pop_front();
        if (m_nodes.isEmpty())
            return;
        m_startLocation = GraphLocation::startOfNode(m_nodes.front());
    }
}


void Path::trimOneBaseOffEnd()
{
    int endPosition = m_endLocation.getPosition();

    //If there is room in the last node to move backward by one, do that.
    if (endPosition > 1)
        m_endLocation.moveLocation(-1);

    //If not, we have to remove the last node and set the Path end to the
    //end of the new last node.
    else
    {
        if (m_nodes.size() > 0)
            m_nodes.pop_back();
        if (m_edges.size() > 0)
            m_edges.pop_back();
        if (m_nodes.isEmpty())
            return;
        m_endLocation = GraphLocation::endOfNode(m_nodes.back());
    }
}



bool Path::containsLocation(GraphLocation location) const
{
    if (location.isNull())
        return false;

    DeBruijnNode * node = location.getNode();

    if (!containsNode(node))
        return false;
    if (containsEntireNode(node))
        return true;

    //If the code got here, then the location's node is partially contained in
    //the path, which means that it is either the first or last node in the
    //path.

    int position = location.getPosition();

    bool isStartNode = node == m_startLocation.getNode();
    bool isEndNode = node == m_endLocation.getNode();

    if (isStartNode && isEndNode)
        return position >= m_startLocation.getPosition() && position <= m_endLocation.getPosition();
    if (isStartNode)
        return position >= m_startLocation.getPosition();
    if (isEndNode)
        return position <= m_endLocation.getPosition();

    return false;
}
