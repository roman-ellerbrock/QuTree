cmake_policy(SET CMP0076 NEW) # enables using cmake V3.12 instead of V3.13

include_directories(${YAML_INCLUDE_DIRS})

add_library(qutreelib)

target_sources(qutreelib PUBLIC
    ${PROJECT_SOURCE_DIR}/src/arithmetics.cpp
    ${PROJECT_SOURCE_DIR}/src/backend/torch_util.cpp
    ${PROJECT_SOURCE_DIR}/src/Graph.cpp
    ${PROJECT_SOURCE_DIR}/src/GraphFactory.cpp
    ${PROJECT_SOURCE_DIR}/src/Operator.cpp
    ${PROJECT_SOURCE_DIR}/src/OperatorFactory.cpp
    ${PROJECT_SOURCE_DIR}/src/TensorNetwork.cpp
)

target_include_directories(qutreelib PUBLIC
    ${PROJECT_SOURCE_DIR}/include
    ${YAML_INCLUDE_DIRS}
    ${NLOHMANN_JSON_DIRS}
)

target_link_libraries(qutreelib PUBLIC
    ${YAML_CPP_LIBRARIES}
    ${TORCH_LIBRARIES}
    nlohmann_json::nlohmann_json
)

set_target_properties(qutreelib 
    PROPERTIES CXX_STANDARD 17
)
