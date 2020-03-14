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


#include "distancedialog.h"
#include "ui_distancedialog.h"

#include "../program/globals.h"
#include "../blast/blastsearch.h"
#include <QMessageBox>
#include "../graph/graphlocation.h"
#include "../graph/assemblygraph.h"
#include "../graph/path.h"
#include "tablewidgetitemint.h"
#include "../program/memory.h"
#include "../program/settings.h"
#include "myprogressdialog.h"
#include "pathsequencecopybutton.h"

DistanceDialog::DistanceDialog(QWidget *parent) :
    QDialog(parent),
    ui(new Ui::DistanceDialog)
{
    ui->setupUi(this);
    setWindowFlags(windowFlags() | Qt::Tool);

    loadSettings();
    searchTypeChanged();
    setInfoTexts();

    ui->resultsTableWidget->setHorizontalHeaderLabels(QStringList() << "Orientation" << "Distance" << "Query 1 path" <<
                                               "Query 2 path" << "Copy distance\nsequence to\nclipboard" << "Distance path");

    //Fill the query combo boxes with any BLAST queries that have paths.
    ui->query1ComboBox->clear();
    ui->query2ComboBox->clear();
    QStringList comboBoxItems;
    for (size_t i = 0; i < g_blastSearch->m_blastQueries.m_queries.size(); ++i)
    {
        if (g_blastSearch->m_blastQueries.m_queries[i]->getPathCount() > 0)
            comboBoxItems.push_back(g_blastSearch->m_blastQueries.m_queries[i]->getName());
    }
    ui->query1ComboBox->addItems(comboBoxItems);
    ui->query2ComboBox->addItems(comboBoxItems);

    //Load the previously used query choices.
    bool rememberedQueriesLoaded = false;
    int indexOfQuery1 = comboBoxItems.indexOf(g_memory->distancePathSearchQuery1);
    int indexOfQuery2 = comboBoxItems.indexOf(g_memory->distancePathSearchQuery2);
    if (indexOfQuery1 != -1 && indexOfQuery2 != -1)
    {
        ui->query1ComboBox->setCurrentIndex(indexOfQuery1);
        ui->query2ComboBox->setCurrentIndex(indexOfQuery2);
        rememberedQueriesLoaded = true;
    }

    //If no queries were successfully loaded, then we try to set the first query
    //to index 0 and the second to index 1, if possible.
    else
    {
        if (ui->query1ComboBox->count() > 0)
            ui->query1ComboBox->setCurrentIndex(0);
        if (ui->query2ComboBox->count() > 1)
            ui->query2ComboBox->setCurrentIndex(1);
    }
    query1Changed();
    query2Changed();

    //If remembered queries were loaded, then load the previously used path
    //choices.
    bool rememberedPathsLoaded = false;
    if (rememberedQueriesLoaded)
    {
        QStringList query1Paths;
        for (int i = 0; i < ui->query1PathComboBox->count(); ++i)
            query1Paths.push_back(ui->query1PathComboBox->itemText(i));
        int indexOfQuery1Path = query1Paths.indexOf(g_memory->distancePathSearchQuery1Path);
        QStringList query2Paths;
        for (int i = 0; i < ui->query2PathComboBox->count(); ++i)
            query2Paths.push_back(ui->query2PathComboBox->itemText(i));
        int indexOfQuery2Path = query2Paths.indexOf(g_memory->distancePathSearchQuery2Path);
        if (indexOfQuery1Path != -1 && indexOfQuery2Path != -1)
        {
            ui->query1PathComboBox->setCurrentIndex(indexOfQuery1Path);
            ui->query2PathComboBox->setCurrentIndex(indexOfQuery2Path);
            rememberedPathsLoaded = true;
        }
    }

    //If the previously used queries and paths were successfully loaded and
    //there are results, display them now.  If not, clear any results that might
    //exist.
    if (rememberedQueriesLoaded && rememberedPathsLoaded &&
            !g_memory->distanceSearchResults.empty())
        fillResultsTable();
    else
        g_memory->distanceSearchResults.clear();

    connect(ui->findPathsButton, SIGNAL(clicked(bool)), this, SLOT(findPaths()));
    connect(ui->query1ComboBox, SIGNAL(currentIndexChanged(int)), this, SLOT(query1Changed()));
    connect(ui->query2ComboBox, SIGNAL(currentIndexChanged(int)), this, SLOT(query2Changed()));
    connect(ui->maxNodesSpinBox, SIGNAL(valueChanged(int)), this, SLOT(saveSettings()));
    connect(ui->maxPathDistanceSpinBox, SIGNAL(valueChanged(int)), this, SLOT(saveSettings()));
    connect(ui->minPathDistanceSpinBox, SIGNAL(valueChanged(int)), this, SLOT(saveSettings()));
    connect(ui->overlapOrientationSameStrandCheckBox, SIGNAL(toggled(bool)), this, SLOT(saveSettings()));
    connect(ui->overlapOrientationOppositeStrandCheckBox, SIGNAL(toggled(bool)), this, SLOT(saveSettings()));
    connect(ui->orientation1CheckBox, SIGNAL(toggled(bool)), this, SLOT(saveSettings()));
    connect(ui->orientation2CheckBox, SIGNAL(toggled(bool)), this, SLOT(saveSettings()));
    connect(ui->orientation3CheckBox, SIGNAL(toggled(bool)), this, SLOT(saveSettings()));
    connect(ui->orientation4CheckBox, SIGNAL(toggled(bool)), this, SLOT(saveSettings()));
    connect(ui->findAllPathsRadioButton, SIGNAL(toggled(bool)), this, SLOT(searchTypeChanged()));
    connect(ui->findAllPathsRadioButton, SIGNAL(toggled(bool)), this, SLOT(saveSettings()));
}

DistanceDialog::~DistanceDialog()
{
    delete ui;
}

void DistanceDialog::findPaths()
{
    if (ui->query1ComboBox->count() < 2)
    {
        QMessageBox::information(this, "Two queries required", "To find paths between queries, you must first conduct a BLAST search "
                                                               "for at least two different queries.");
        return;
    }

    if (!g_settings->distanceOverlapOrientationSameStrand &&
            !g_settings->distanceOverlapOrientationOppositeStrand &&
            !g_settings->distanceOrientation1 &&
            !g_settings->distanceOrientation2 &&
            !g_settings->distanceOrientation3 &&
            !g_settings->distanceOrientation4)
    {
        QMessageBox::information(this, "Orientation required", "At least one of the six possible orientations must be ticked "
                                                               "before finding paths between queries.");
        return;
    }

    BlastQuery * query1 = g_blastSearch->m_blastQueries.getQueryFromName(ui->query1ComboBox->currentText());
    BlastQuery * query2 = g_blastSearch->m_blastQueries.getQueryFromName(ui->query2ComboBox->currentText());

    if (query1 == 0 || query2 == 0)
        return;

    if (query1 == query2)
    {
        QMessageBox::information(this, "Same query", "The two selected queries are the same. To find paths between queries, you must select two different queries.");
        return;
    }

    //Remember which queries and paths were used for this search.
    g_memory->distancePathSearchQuery1 = ui->query1ComboBox->currentText();
    g_memory->distancePathSearchQuery2 = ui->query2ComboBox->currentText();
    g_memory->distancePathSearchQuery1Path = ui->query1PathComboBox->currentText();
    g_memory->distancePathSearchQuery2Path = ui->query2PathComboBox->currentText();

    //If the current index of the path combo box is 0, that means either there is
    //only one path and it is selected or there are multiple paths and "all" is
    //selected.  In either case, we can use all paths.
    QList<BlastQueryPath> query1Paths;
    if (ui->query1PathComboBox->currentIndex() == 0)
        query1Paths = query1->getPaths();

    //If the current index is not 0, then there are multiple paths and one path
    //in particular is selected.
    else
        query1Paths.push_back(query1->getPaths()[ui->query1PathComboBox->currentIndex() - 1]);

    //Repeat for query 2.
    QList<BlastQueryPath> query2Paths;
    if (ui->query2PathComboBox->currentIndex() == 0)
        query2Paths = query2->getPaths();
    else
        query2Paths.push_back(query2->getPaths()[ui->query2PathComboBox->currentIndex() - 1]);

    //If an 'all paths search is conducted, then we must make sure that the min
    //distance is less than the max distance.
    if (g_settings->findAllDistancePaths)
    {
        if (g_settings->minDistancePathLength > g_settings->maxDistancePathLength)
        {
            QMessageBox::information(this, "Min greater than max", "The minimum path distance must be less than or equal to the maximum path distance.");
            return;
        }
    }

    //Run the path search.  This is in a separate code block so the progress
    //dialog is destroyed when the search is finished.
    {
        MyProgressDialog progress(this, "Finding paths between queries...", false);
        progress.setWindowModality(Qt::WindowModal);
        progress.show();

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
    }

    fillResultsTable();

    if (g_memory->distanceSearchResults.size() == 0)
    {
        if (g_settings->findAllDistancePaths)
            QMessageBox::information(this, "No paths", "No paths were found between the two given queries with the selected orientations and settings.");
        else
            QMessageBox::information(this, "No path", "No path exists between the two given queries with the selected orientations.");
    }
}


//These functions populate the path combo boxes for each query.
void DistanceDialog::query1Changed()
{
    BlastQuery * query1 = g_blastSearch->m_blastQueries.getQueryFromName(ui->query1ComboBox->currentText());
    if (query1 == 0)
        return;

    fillPathComboBox(query1, ui->query1PathComboBox);
}
void DistanceDialog::query2Changed()
{
    BlastQuery * query2 = g_blastSearch->m_blastQueries.getQueryFromName(ui->query2ComboBox->currentText());
    if (query2 == 0)
        return;

    fillPathComboBox(query2, ui->query2PathComboBox);
}
void DistanceDialog::fillPathComboBox(BlastQuery * query, QComboBox * comboBox)
{
    comboBox->clear();

    QStringList comboBoxItems;
    QList<BlastQueryPath> paths = query->getPaths();

    if (paths.size() == 0)
        return;
    else if (paths.size() == 1)
        comboBoxItems.push_back(paths[0].getPath().getString(true));
    else
    {
        comboBoxItems.push_back("all");
        for (int i = 0; i < paths.size(); ++i)
            comboBoxItems.push_back(paths[i].getPath().getString(true));
    }

    comboBox->addItems(comboBoxItems);
}



void DistanceDialog::fillResultsTable()
{
    int pathCount = g_memory->distanceSearchResults.size();

    ui->resultsTableWidget->clearContents();
    ui->resultsTableWidget->setSortingEnabled(false);
    ui->resultsTableWidget->setRowCount(pathCount);

    for (int i = 0; i < pathCount; ++i)
    {
        QTableWidgetItem * query1Path = new QTableWidgetItem(g_memory->distanceSearchResults[i].m_query1Path.getString(true));
        query1Path->setFlags(Qt::ItemIsEnabled | Qt::ItemIsSelectable);

        QTableWidgetItem * query2Path = new QTableWidgetItem(g_memory->distanceSearchResults[i].m_query2Path.getString(true));
        query2Path->setFlags(Qt::ItemIsEnabled | Qt::ItemIsSelectable);

        QTableWidgetItem * orientation = new QTableWidgetItem(g_memory->distanceSearchResults[i].m_orientation);
        orientation->setFlags(Qt::ItemIsEnabled | Qt::ItemIsSelectable);

        int distanceValue = g_memory->distanceSearchResults[i].m_distance;
        TableWidgetItemInt * distance = new TableWidgetItemInt(formatIntForDisplay(distanceValue), distanceValue);
        distance->setFlags(Qt::ItemIsEnabled | Qt::ItemIsSelectable);

        QTableWidgetItem * distancePath = new QTableWidgetItem(g_memory->distanceSearchResults[i].m_distancePath.getString(true));
        distancePath->setFlags(Qt::ItemIsEnabled | Qt::ItemIsSelectable);

        QByteArray pathSequence = g_memory->distanceSearchResults[i].m_distancePath.getPathSequence();
        QString pathStart;
        if (pathSequence.length() <= 8)
            pathStart = pathSequence;
        else
            pathStart = pathSequence.left(8) + "...";

        QTableWidgetItem * sequenceCopy = new QTableWidgetItem(pathStart);
        PathSequenceCopyButton * sequenceCopyButton = new PathSequenceCopyButton(pathSequence, pathStart);

        ui->resultsTableWidget->setItem(i, 0, orientation);
        ui->resultsTableWidget->setItem(i, 1, distance);
        ui->resultsTableWidget->setItem(i, 2, query1Path);
        ui->resultsTableWidget->setItem(i, 3, query2Path);
        ui->resultsTableWidget->setItem(i, 4, sequenceCopy);
        ui->resultsTableWidget->setCellWidget(i, 4, sequenceCopyButton);
        ui->resultsTableWidget->setItem(i, 5, distancePath);
    }

    ui->resultsTableWidget->resizeColumns();
    ui->resultsTableWidget->setSortingEnabled(true);
}



void DistanceDialog::loadSettings()
{
    ui->maxNodesSpinBox->setValue(g_settings->distancePathSearchDepth + 1);
    ui->maxNodesSpinBox->setMinimum(g_settings->distancePathSearchDepth.min);
    ui->maxNodesSpinBox->setMaximum(g_settings->distancePathSearchDepth.max);

    ui->minPathDistanceSpinBox->setValue(g_settings->minDistancePathLength);
    ui->minPathDistanceSpinBox->setMinimum(g_settings->minDistancePathLength.min);
    ui->minPathDistanceSpinBox->setMaximum(g_settings->minDistancePathLength.max);

    ui->maxPathDistanceSpinBox->setValue(g_settings->maxDistancePathLength);
    ui->maxPathDistanceSpinBox->setMinimum(g_settings->maxDistancePathLength.min);
    ui->maxPathDistanceSpinBox->setMaximum(g_settings->maxDistancePathLength.max);

    ui->overlapOrientationSameStrandCheckBox->setChecked(g_settings->distanceOverlapOrientationSameStrand);
    ui->overlapOrientationOppositeStrandCheckBox->setChecked(g_settings->distanceOverlapOrientationOppositeStrand);
    ui->orientation1CheckBox->setChecked(g_settings->distanceOrientation1);
    ui->orientation2CheckBox->setChecked(g_settings->distanceOrientation2);
    ui->orientation3CheckBox->setChecked(g_settings->distanceOrientation3);
    ui->orientation4CheckBox->setChecked(g_settings->distanceOrientation4);
    ui->findAllPathsRadioButton->setChecked(g_settings->findAllDistancePaths);
    ui->findShortestPathsRadioButton->setChecked(!g_settings->findAllDistancePaths);
}

void DistanceDialog::saveSettings()
{
    g_settings->distancePathSearchDepth = ui->maxNodesSpinBox->value() - 1;
    g_settings->minDistancePathLength = ui->minPathDistanceSpinBox->value();
    g_settings->maxDistancePathLength = ui->maxPathDistanceSpinBox->value();
    g_settings->distanceOverlapOrientationSameStrand = ui->overlapOrientationSameStrandCheckBox->isChecked();
    g_settings->distanceOverlapOrientationOppositeStrand = ui->overlapOrientationOppositeStrandCheckBox->isChecked();
    g_settings->distanceOrientation1 = ui->orientation1CheckBox->isChecked();
    g_settings->distanceOrientation2 = ui->orientation2CheckBox->isChecked();
    g_settings->distanceOrientation3 = ui->orientation3CheckBox->isChecked();
    g_settings->distanceOrientation4 = ui->orientation4CheckBox->isChecked();
    g_settings->findAllDistancePaths = ui->findAllPathsRadioButton->isChecked();
}


void DistanceDialog::setInfoTexts()
{
    ui->queriesInfoText->setInfoText("Select two different BLAST queries here for which you want to find "
                                     "the distance between. Queries are only available if Bandage has "
                                     "found at least one graph path.<br><br>"
                                     "If a query has more than one graph path, then you can use "
                                     "all of the paths in the search or you can select only one.");

    ui->maxNodesInfoText->setInfoText("This is the maximum number of nodes that can be in a path "
                                      "between the two queries.<br><br>"
                                      "Larger values can allow the search to find more complex "
                                      "paths in the graph but at a performance cost.");

    ui->minPathDistanceInfoText->setInfoText("Paths shorter than this length (measured in base "
                                             "pairs) will not be included in the search results.");

    ui->maxPathDistanceInfoText->setInfoText("Paths longer than this length (measured in base "
                                             "pairs) will not be included in the search results.");

    ui->orientationsInfoText->setInfoText("Each possible query orientation can be included or "
                                          "excluded from the search using these tick boxes:"
                                          "<ul>"
                                          "<li>1-&#62; 2-&#62; The two queries are on the same strand "
                                          "of DNA, with query 1 occurring upstream of query 2.</li>"
                                          "<li>2-&#62; 1-&#62; The two queries are on the same strand "
                                          "of DNA, with query 1 occurring downstream of query 2.</li>"
                                          "<li>1-&#62; &#60;-2 The two queries are on different strands "
                                          "of DNA, with their 3' ends closer than their 5' ends.</li>"
                                          "<li>&#60;-1 2-&#62; The two queries are on different strands "
                                          "of DNA, with their 5' ends closer than their 3' ends.</li>"
                                          "</ul>");

    ui->findShortestPathsInfoText->setInfoText("If this option is chosen, then the path search will only "
                                               "look for one shortest path between two queries. If no "
                                               "paths are found, then the two queries are not connected "
                                               "in the graph.<br><br>"
                                               "The returned path is guaranteed to be the shortest "
                                               "path, though it is possible for other equally short paths "
                                               "to be present in the graph.<br><br>"
                                               "This search examines the entire graph and the settings "
                                               "to the right are not used.");
    ui->findAllPathsInfoText->setInfoText("If this option is chosen, then the path search will look for "
                                          "all possible paths between the two queries.<br><br>"
                                          "Since there could be many (possible infinite) paths between "
                                          "two queries, the search is limited by using the settings to "
                                          "the right.");
}



void DistanceDialog::searchTypeChanged()
{
    ui->maxNodesInfoText->setEnabled(ui->findAllPathsRadioButton->isChecked());
    ui->maxNodesLabel->setEnabled(ui->findAllPathsRadioButton->isChecked());
    ui->maxNodesSpinBox->setEnabled(ui->findAllPathsRadioButton->isChecked());
    ui->minPathDistanceInfoText->setEnabled(ui->findAllPathsRadioButton->isChecked());
    ui->minPathDistanceLabel->setEnabled(ui->findAllPathsRadioButton->isChecked());
    ui->minPathDistanceSpinBox->setEnabled(ui->findAllPathsRadioButton->isChecked());
    ui->maxPathDistanceInfoText->setEnabled(ui->findAllPathsRadioButton->isChecked());
    ui->maxPathDistanceLabel->setEnabled(ui->findAllPathsRadioButton->isChecked());
    ui->maxPathDistanceSpinBox->setEnabled(ui->findAllPathsRadioButton->isChecked());
}
