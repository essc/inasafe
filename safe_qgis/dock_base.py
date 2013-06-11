# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file 'dock_base.ui'
#
# Created: Tue Jun 11 10:17:52 2013
#      by: PyQt4 UI code generator 4.10
#
# WARNING! All changes made in this file will be lost!

from PyQt4 import QtCore, QtGui

try:
    _fromUtf8 = QtCore.QString.fromUtf8
except AttributeError:
    def _fromUtf8(s):
        return s

try:
    _encoding = QtGui.QApplication.UnicodeUTF8
    def _translate(context, text, disambig):
        return QtGui.QApplication.translate(context, text, disambig, _encoding)
except AttributeError:
    def _translate(context, text, disambig):
        return QtGui.QApplication.translate(context, text, disambig)

class Ui_DockBase(object):
    def setupUi(self, DockBase):
        DockBase.setObjectName(_fromUtf8("DockBase"))
        DockBase.resize(393, 547)
        icon = QtGui.QIcon()
        icon.addPixmap(QtGui.QPixmap(_fromUtf8(":/plugins/inasafe/icon.png")), QtGui.QIcon.Normal, QtGui.QIcon.Off)
        DockBase.setWindowIcon(icon)
        self.dockWidgetContents = QtGui.QWidget()
        self.dockWidgetContents.setObjectName(_fromUtf8("dockWidgetContents"))
        self.gridLayout = QtGui.QGridLayout(self.dockWidgetContents)
        self.gridLayout.setContentsMargins(3, 0, 3, 3)
        self.gridLayout.setObjectName(_fromUtf8("gridLayout"))
        self.verticalLayout_4 = QtGui.QVBoxLayout()
        self.verticalLayout_4.setContentsMargins(3, -1, -1, -1)
        self.verticalLayout_4.setObjectName(_fromUtf8("verticalLayout_4"))
        self.wvResults = MessageViewer(self.dockWidgetContents)
        sizePolicy = QtGui.QSizePolicy(QtGui.QSizePolicy.Preferred, QtGui.QSizePolicy.MinimumExpanding)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.wvResults.sizePolicy().hasHeightForWidth())
        self.wvResults.setSizePolicy(sizePolicy)
        self.wvResults.setMinimumSize(QtCore.QSize(0, 50))
        self.wvResults.setProperty("url", QtCore.QUrl(_fromUtf8("about:blank")))
        self.wvResults.setObjectName(_fromUtf8("wvResults"))
        self.verticalLayout_4.addWidget(self.wvResults)
        self.horizontalLayout_3 = QtGui.QHBoxLayout()
        self.horizontalLayout_3.setObjectName(_fromUtf8("horizontalLayout_3"))
        self.label_9 = QtGui.QLabel(self.dockWidgetContents)
        self.label_9.setMinimumSize(QtCore.QSize(64, 64))
        self.label_9.setMaximumSize(QtCore.QSize(64, 64))
        self.label_9.setText(_fromUtf8(""))
        self.label_9.setPixmap(QtGui.QPixmap(_fromUtf8(":/plugins/inasafe/bnpb_logo_64.png")))
        self.label_9.setScaledContents(True)
        self.label_9.setAlignment(QtCore.Qt.AlignCenter)
        self.label_9.setObjectName(_fromUtf8("label_9"))
        self.horizontalLayout_3.addWidget(self.label_9)
        self.verticalLayout = QtGui.QVBoxLayout()
        self.verticalLayout.setObjectName(_fromUtf8("verticalLayout"))
        self.label_3 = QtGui.QLabel(self.dockWidgetContents)
        font = QtGui.QFont()
        font.setBold(True)
        font.setWeight(75)
        self.label_3.setFont(font)
        self.label_3.setAlignment(QtCore.Qt.AlignCenter)
        self.label_3.setObjectName(_fromUtf8("label_3"))
        self.verticalLayout.addWidget(self.label_3)
        self.horizontalLayout_2 = QtGui.QHBoxLayout()
        self.horizontalLayout_2.setObjectName(_fromUtf8("horizontalLayout_2"))
        spacerItem = QtGui.QSpacerItem(40, 20, QtGui.QSizePolicy.Expanding, QtGui.QSizePolicy.Minimum)
        self.horizontalLayout_2.addItem(spacerItem)
        self.label_4 = QtGui.QLabel(self.dockWidgetContents)
        self.label_4.setText(_fromUtf8(""))
        self.label_4.setPixmap(QtGui.QPixmap(_fromUtf8(":/plugins/inasafe/Australian-AID-Identifier_blue-red-small.png")))
        self.label_4.setAlignment(QtCore.Qt.AlignCenter)
        self.label_4.setWordWrap(True)
        self.label_4.setObjectName(_fromUtf8("label_4"))
        self.horizontalLayout_2.addWidget(self.label_4)
        spacerItem1 = QtGui.QSpacerItem(40, 20, QtGui.QSizePolicy.Expanding, QtGui.QSizePolicy.Minimum)
        self.horizontalLayout_2.addItem(spacerItem1)
        self.label = QtGui.QLabel(self.dockWidgetContents)
        self.label.setText(_fromUtf8(""))
        self.label.setPixmap(QtGui.QPixmap(_fromUtf8(":/plugins/inasafe/wb-logo-small.png")))
        self.label.setObjectName(_fromUtf8("label"))
        self.horizontalLayout_2.addWidget(self.label)
        self.verticalLayout.addLayout(self.horizontalLayout_2)
        self.horizontalLayout_3.addLayout(self.verticalLayout)
        self.verticalLayout_4.addLayout(self.horizontalLayout_3)
        self.gridLayout.addLayout(self.verticalLayout_4, 2, 0, 1, 1)
        self.grpQuestion = QtGui.QGroupBox(self.dockWidgetContents)
        self.grpQuestion.setObjectName(_fromUtf8("grpQuestion"))
        self.gridLayout_3 = QtGui.QGridLayout(self.grpQuestion)
        self.gridLayout_3.setContentsMargins(0, 6, 0, 0)
        self.gridLayout_3.setVerticalSpacing(1)
        self.gridLayout_3.setObjectName(_fromUtf8("gridLayout_3"))
        self.cboHazard = QtGui.QComboBox(self.grpQuestion)
        self.cboHazard.setInsertPolicy(QtGui.QComboBox.InsertAlphabetically)
        self.cboHazard.setObjectName(_fromUtf8("cboHazard"))
        self.gridLayout_3.addWidget(self.cboHazard, 0, 0, 1, 2)
        self.label_7 = QtGui.QLabel(self.grpQuestion)
        self.label_7.setObjectName(_fromUtf8("label_7"))
        self.gridLayout_3.addWidget(self.label_7, 1, 0, 1, 1)
        self.cboExposure = QtGui.QComboBox(self.grpQuestion)
        self.cboExposure.setInsertPolicy(QtGui.QComboBox.InsertAlphabetically)
        self.cboExposure.setObjectName(_fromUtf8("cboExposure"))
        self.gridLayout_3.addWidget(self.cboExposure, 3, 0, 1, 2)
        self.label_8 = QtGui.QLabel(self.grpQuestion)
        self.label_8.setObjectName(_fromUtf8("label_8"))
        self.gridLayout_3.addWidget(self.label_8, 4, 0, 1, 1)
        self.cboFunction = QtGui.QComboBox(self.grpQuestion)
        sizePolicy = QtGui.QSizePolicy(QtGui.QSizePolicy.Preferred, QtGui.QSizePolicy.Fixed)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.cboFunction.sizePolicy().hasHeightForWidth())
        self.cboFunction.setSizePolicy(sizePolicy)
        self.cboFunction.setInsertPolicy(QtGui.QComboBox.InsertAlphabetically)
        self.cboFunction.setObjectName(_fromUtf8("cboFunction"))
        self.gridLayout_3.addWidget(self.cboFunction, 5, 0, 1, 1)
        self.cboAggregation = QtGui.QComboBox(self.grpQuestion)
        self.cboAggregation.setObjectName(_fromUtf8("cboAggregation"))
        self.gridLayout_3.addWidget(self.cboAggregation, 7, 0, 1, 2)
        self.toolFunctionOptions = QtGui.QToolButton(self.grpQuestion)
        self.toolFunctionOptions.setEnabled(False)
        sizePolicy = QtGui.QSizePolicy(QtGui.QSizePolicy.Fixed, QtGui.QSizePolicy.Fixed)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.toolFunctionOptions.sizePolicy().hasHeightForWidth())
        self.toolFunctionOptions.setSizePolicy(sizePolicy)
        self.toolFunctionOptions.setMinimumSize(QtCore.QSize(27, 27))
        self.toolFunctionOptions.setMaximumSize(QtCore.QSize(16777215, 27))
        icon1 = QtGui.QIcon()
        icon1.addPixmap(QtGui.QPixmap(_fromUtf8(":/plugins/inasafe/edit.png")), QtGui.QIcon.Normal, QtGui.QIcon.On)
        self.toolFunctionOptions.setIcon(icon1)
        self.toolFunctionOptions.setIconSize(QtCore.QSize(27, 27))
        self.toolFunctionOptions.setObjectName(_fromUtf8("toolFunctionOptions"))
        self.gridLayout_3.addWidget(self.toolFunctionOptions, 5, 1, 1, 1)
        self.label_2 = QtGui.QLabel(self.grpQuestion)
        self.label_2.setObjectName(_fromUtf8("label_2"))
        self.gridLayout_3.addWidget(self.label_2, 6, 0, 1, 1)
        self.gridLayout.addWidget(self.grpQuestion, 0, 0, 1, 1)
        self.horizontalLayout = QtGui.QHBoxLayout()
        self.horizontalLayout.setObjectName(_fromUtf8("horizontalLayout"))
        self.pbnHelp = QtGui.QPushButton(self.dockWidgetContents)
        self.pbnHelp.setObjectName(_fromUtf8("pbnHelp"))
        self.horizontalLayout.addWidget(self.pbnHelp)
        spacerItem2 = QtGui.QSpacerItem(40, 20, QtGui.QSizePolicy.Expanding, QtGui.QSizePolicy.Minimum)
        self.horizontalLayout.addItem(spacerItem2)
        self.pbnPrint = QtGui.QPushButton(self.dockWidgetContents)
        self.pbnPrint.setObjectName(_fromUtf8("pbnPrint"))
        self.horizontalLayout.addWidget(self.pbnPrint)
        spacerItem3 = QtGui.QSpacerItem(40, 20, QtGui.QSizePolicy.Expanding, QtGui.QSizePolicy.Minimum)
        self.horizontalLayout.addItem(spacerItem3)
        self.pbnRunStop = QtGui.QPushButton(self.dockWidgetContents)
        self.pbnRunStop.setObjectName(_fromUtf8("pbnRunStop"))
        self.horizontalLayout.addWidget(self.pbnRunStop)
        self.gridLayout.addLayout(self.horizontalLayout, 5, 0, 1, 1)
        DockBase.setWidget(self.dockWidgetContents)
        self.label_7.setBuddy(self.cboExposure)
        self.label_8.setBuddy(self.cboFunction)

        self.retranslateUi(DockBase)
        QtCore.QMetaObject.connectSlotsByName(DockBase)
        DockBase.setTabOrder(self.cboHazard, self.cboExposure)
        DockBase.setTabOrder(self.cboExposure, self.cboFunction)
        DockBase.setTabOrder(self.cboFunction, self.pbnRunStop)
        DockBase.setTabOrder(self.pbnRunStop, self.pbnHelp)

    def retranslateUi(self, DockBase):
        DockBase.setWindowTitle(_translate("DockBase", "InaSAFE", None))
        self.label_3.setText(_translate("DockBase", "Supported by:", None))
        self.grpQuestion.setTitle(_translate("DockBase", "Question: In the event of", None))
        self.label_7.setText(_translate("DockBase", "How many", None))
        self.label_8.setText(_translate("DockBase", "&Might", None))
        self.toolFunctionOptions.setToolTip(_translate("DockBase", "Configure Impact Function Parameter", None))
        self.toolFunctionOptions.setText(_translate("DockBase", "...", None))
        self.label_2.setText(_translate("DockBase", "Aggregate results by", None))
        self.pbnHelp.setText(_translate("DockBase", "Help", None))
        self.pbnPrint.setText(_translate("DockBase", "Print...", None))
        self.pbnRunStop.setText(_translate("DockBase", "Run", None))

from message_viewer import MessageViewer
import resources_rc
