#include "mainwindow.h"
#include <QtCore/QDir>
#include <QtCore/QFileInfo>
#include <QtCore/QByteArray>
#include <QtCore/QStandardPaths>
#include <QtCore/QRect>
#include <QtGui/QScreen>
#include <QtWidgets/QApplication>
#include <QtWidgets/QStyleFactory>
#include <QtWidgets/QMessageBox>
#include <QtCore/QTimer>

// This is where our application starts - the main entry point
// We set up the Qt application, configure the look and feel, and launch the
// main window
int main(int argc, char *argv[]) {
  // Note: Avoid forcing Qt plugin/runtime paths here to prevent mixing
  // different Qt installations. Use deployment (windeployqt) or environment
  // variables outside the app to control plugin loading.
  // Create the Qt application object - this manages the entire GUI application
  QApplication app(argc, argv);

  // Set up basic information about our application
  // This helps with things like window titles, taskbar icons, and system
  // integration
  app.setApplicationName("Kernel Function Comparison Tool");
  app.setApplicationVersion("1.0.0");
  app.setOrganizationName("PDIP Research");
  app.setApplicationDisplayName(
      "Primal-Dual Interior Point Methods Comparison");

  // Choose Fusion for consistency and load QSS theme from resources
  app.setStyle(QStyleFactory::create("Fusion"));
  QFile themeFile(":/qss/qss/app_theme.qss");
  if (themeFile.open(QIODevice::ReadOnly | QIODevice::Text)) {
    const QString qss = QString::fromUtf8(themeFile.readAll());
    qApp->setStyleSheet(qss);
  }

  // Create the main window of our application
  // This is where all the magic happens - the user interface for comparing
  // kernel functions
  MainWindow window;

  // Show window
  window.show();

  // Safety: if window is off-screen, move it to a safe location
  QRect screenGeom = QApplication::primaryScreen()->geometry();
  QRect winGeom = window.frameGeometry();
  if (!screenGeom.intersects(winGeom)) {
    window.move(100, 100);
  }

  // Start the application's event loop
  // This keeps the application running and responsive to user interactions
  // The app will keep running until the user closes the main window
  return app.exec();
}