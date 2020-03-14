#include "pathsequencecopybutton.h"

#include <QApplication>
#include <QClipboard>

PathSequenceCopyButton::PathSequenceCopyButton(QByteArray pathSequence, QString pathStart) :
    m_pathSequence(pathSequence)
{
    setText(pathStart);
    connect(this, SIGNAL(clicked(bool)), this, SLOT(copySequenceToClipboard()));
}


void PathSequenceCopyButton::copySequenceToClipboard()
{
    QClipboard * clipboard = QApplication::clipboard();
    clipboard->setText(m_pathSequence);
}
