# CMake configuration for the Kernel Function Optimization Platform
# This file tells the build system how to compile our C++ application

# We need CMake version 3.16 or higher for Qt6 support
cmake_minimum_required(VERSION 3.16)

# Define our project - this is the main application for comparing kernel functions
project(KernelFunction VERSION 1.0.0 LANGUAGES CXX)

# Set the C++ standard to C++17
# This gives us access to modern C++ features like std::optional, std::variant, etc.
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Find Qt6 components that we need
# Core: Basic Qt functionality
# Widgets: GUI components (buttons, windows, etc.)
# Charts: For displaying performance charts (optional, can be added later)
find_package(Qt6 REQUIRED COMPONENTS Core Widgets Concurrent)
find_package(Qt6 COMPONENTS Charts QUIET)
find_package(Eigen3 REQUIRED)

# Enable Qt's automatic code generation
# MOC: Meta-Object Compiler - handles Qt's signal/slot system
# UIC: User Interface Compiler - processes .ui files
# RCC: Resource Compiler - handles embedded resources
set(CMAKE_AUTOMOC ON)
set(CMAKE_AUTOUIC ON)
set(CMAKE_AUTORCC ON)

# Tell the compiler where to find our header files
# This allows us to use #include "filename.h" instead of #include "src/filename.h"
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/src)

# === GUI APPLICATION SOURCE FILES ===
# These are the files needed to build the graphical user interface version

set(GUI_SOURCES
    src/main.cpp                    # Application entry point
    src/mainwindow.cpp              # Main window implementation
    src/gui/algorithm_dialog.cpp    # Algorithm details dialog
    src/gui/matrix_input_widget.cpp # Matrix input widget
    src/ui/kernel_manager.cpp       # Data-driven kernel registry
    src/ui/parameter_editor_widget.cpp # Dynamic parameter editor
    src/engine/ipm_engine.cpp       # IPM engine
    src/algorithms/algorithm_base.cpp    # Base algorithm class
    src/algorithms/algorithm2.cpp        # Algorithm 2 implementation
    src/kernels/kernel_base.cpp          # Base kernel class
    src/kernels/trigonometric_kernel.cpp      # Trigonometric kernel
    src/kernels/exponential_parametric_kernel.cpp  # Exponential-parametric kernel
    src/kernels/parameterized_log_kernel.cpp       # Parameterized log kernel
    src/kernels/parametric_family_kernel.cpp       # Parametric family kernel
    src/kernels/log_exponential_kernel.cpp         # Log-exponential kernel
    src/kernels/bachir_kernel.cpp                  # Bachir kernel
    src/kernels/touil_hyperbolic_kernel.cpp        # Touil22 kernel
    src/kernels/ijnna23_fractional_kernel.cpp      # IJNAA23 kernel
    src/kernels/wu_hyperbolic_kernel.cpp           # Wu25 kernel
    src/kernels/bachir_convex_combo_kernel.cpp     # Bachir25 convex combo kernel
    src/kernels/derbal20_param_log_kernel.cpp      # Derbal20 param log kernel
    src/utils/test_data.cpp                        # Test problem data
    resources/app_resources.qrc                    # Embedded resources (QSS)
)

# Optionally add charts source if available
if (Qt6Charts_FOUND)
  list(APPEND GUI_SOURCES src/gui/cpu_time_chart_widget.cpp)
  set(HAVE_QT_CHARTS TRUE)
endif()

# Header files for the GUI application
set(GUI_HEADERS
    src/mainwindow.h                # Main window class definition
    src/gui/algorithm_dialog.h      # Algorithm details dialog header
    src/gui/matrix_input_widget.h   # Matrix input widget header
    src/ui/kernel_manager.h         # Kernel registry
    src/ui/parameter_editor_widget.h # Parameter editor widget
    src/engine/ipm_engine.h         # IPM engine
    src/algorithms/algorithm_base.h      # Base algorithm class header
    src/algorithms/algorithm2.h          # Algorithm 2 header
    src/kernels/kernel_base.h            # Base kernel class header
    src/kernels/trigonometric_kernel.h        # Trigonometric kernel header
    src/kernels/exponential_parametric_kernel.h  # Exponential-parametric kernel header
    src/kernels/parameterized_log_kernel.h       # Parameterized log kernel header
    src/kernels/parametric_family_kernel.h       # Parametric family kernel header
    src/kernels/log_exponential_kernel.h         # Log-exponential kernel header
    src/kernels/bachir_kernel.h                  # Bachir kernel header
    src/kernels/touil_hyperbolic_kernel.h        # Touil22 kernel header
    src/kernels/ijnna23_fractional_kernel.h      # IJNAA23 kernel header
    src/kernels/wu_hyperbolic_kernel.h           # Wu25 kernel header
    src/kernels/bachir_convex_combo_kernel.h     # Bachir25 convex combo kernel header
    src/kernels/derbal20_param_log_kernel.h      # Derbal20 param log kernel header
    src/utils/test_data.h                        # Test data header
)

# === CONSOLE APPLICATION SOURCE FILES ===
# These are the files needed to build the command-line version
# Useful for testing and batch processing

set(CONSOLE_SOURCES
    src/main_console.cpp            # Console application entry point
    src/engine/ipm_engine.cpp       # IPM engine
    src/algorithms/algorithm_base.cpp    # Base algorithm class
    src/algorithms/algorithm1.cpp        # Algorithm 1 implementation
    src/algorithms/algorithm2.cpp        # Algorithm 2 implementation
    src/algorithms/algorithm3.cpp        # Algorithm 3 implementation
    src/algorithms/algorithm4.cpp        # Algorithm 4 implementation
    src/algorithms/algorithm5.cpp        # Algorithm 5 implementation
    src/kernels/kernel_base.cpp          # Base kernel class
    src/kernels/trigonometric_kernel.cpp      # Trigonometric kernel
    src/kernels/exponential_parametric_kernel.cpp  # Exponential-parametric kernel
    src/kernels/parameterized_log_kernel.cpp       # Parameterized log kernel
    src/kernels/parametric_family_kernel.cpp       # Parametric family kernel
    src/kernels/log_exponential_kernel.cpp         # Log-exponential kernel
    src/kernels/bachir_kernel.cpp                  # Bachir kernel
    src/kernels/touil_hyperbolic_kernel.cpp        # Touil22 kernel
    src/kernels/ijnna23_fractional_kernel.cpp      # IJNAA23 kernel
    src/kernels/wu_hyperbolic_kernel.cpp           # Wu25 kernel
    src/kernels/bachir_convex_combo_kernel.cpp     # Bachir25 convex combo kernel
    src/kernels/derbal20_param_log_kernel.cpp      # Derbal20 param log kernel
    src/utils/test_data.cpp                        # Test problem data
)

# === CREATE THE EXECUTABLES ===

# Build the GUI application - this creates the main program with windows and buttons
add_executable(kernel_ui ${GUI_SOURCES} ${GUI_HEADERS})

# Build the console application - this creates a command-line version
add_executable(kernel_console ${CONSOLE_SOURCES})

# === LINK LIBRARIES ===

# Link Qt libraries to the GUI application
# These provide the GUI functionality (windows, buttons, charts, etc.)
target_link_libraries(kernel_ui 
    Qt6::Core      # Basic Qt functionality
    Qt6::Widgets   # GUI components
    Qt6::Concurrent # QtConcurrent for background tasks
    Eigen3::Eigen
)

# Link Qt Charts only if found (graceful fallback)
if (Qt6Charts_FOUND)
  target_compile_definitions(kernel_ui PRIVATE CPU_CHARTS_ENABLED)
  target_link_libraries(kernel_ui Qt6::Charts)
  message(STATUS "Qt Charts found: CPU Time graphs enabled")
else()
  message(WARNING "Qt Charts NOT found: CPU Time graphs will be disabled. Install Qt Charts add-on for Qt 6.7.3.")
endif()

# Link Eigen to console application for headers and transitive includes
target_link_libraries(kernel_console 
    Eigen3::Eigen
)

# === COMPILER OPTIONS ===
# Set warning levels and other compiler-specific options

if(MSVC)
    # Microsoft Visual C++ compiler options
    # /W4 enables most warnings to help catch potential issues
    target_compile_options(kernel_ui PRIVATE /W4)
    target_compile_options(kernel_console PRIVATE /W4)
else()
    # GCC/Clang compiler options
    # -Wall: Enable most warnings
    # -Wextra: Enable extra warnings
    # -Wpedantic: Strict ISO C++ compliance warnings
    target_compile_options(kernel_ui PRIVATE -Wall -Wextra -Wpedantic)
    target_compile_options(kernel_console PRIVATE -Wall -Wextra -Wpedantic)
endif()

# === WINDOWS-SPECIFIC SETTINGS ===

if(WIN32)
    # Set properties for the GUI application on Windows
    set_target_properties(kernel_ui PROPERTIES
        WIN32_EXECUTABLE FALSE  # Don't hide the console window (useful for debugging)
    )
    
    # Fix linking issues with MinGW (GNU compiler on Windows)
    if(MINGW)
        target_link_libraries(kernel_ui mingw32)
    endif()
endif()

# === INSTALLATION ===
# Define where the compiled programs should be installed
# This is used when someone runs "make install" or "cmake --install"

install(TARGETS kernel_ui kernel_console DESTINATION bin) 