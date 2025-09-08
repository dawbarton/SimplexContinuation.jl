# Test Organization for SimplexContinuation.jl

This directory contains the comprehensive test suite for the SimplexContinuation.jl package. The tests are organized into separate files for better maintainability and clarity.

## Test Structure

### Main Test Runner
- **`runtests.jl`** - Main test runner that includes all unit test files

### Unit Test Files (`unit/` directory)

#### Core Functionality Tests
- **`test_simplex_constructors.jl`** - Tests for all Simplex constructor variants
  - Basic construction from vertices
  - Type-specific construction
  - Uninitialized construction
  - Copy constructors and type conversion
  - Error handling for invalid inputs

- **`test_simplex_properties.jl`** - Tests for simplex property and utility functions
  - `simplex_dimension()` and `space_dimension()`
  - `eltype()` function
  - `sort!()` function with various data types
  - Type consistency across operations
  - Memory safety verification

#### Geometric Operations Tests
- **`test_geometric_reflection.jl`** - Tests for general geometric reflection
  - 2D and 3D geometric reflections
  - Reflection mathematical properties (involution)
  - Error handling (invalid facet indices, non-full-dimensional)
  - Different numeric types (rational, BigFloat, etc.)
  - Geometric correctness verification

#### Freudenthal-Specific Tests
- **`test_freudenthal.jl`** - Tests for Freudenthal triangulation functions
  - `freudenthal_initial_simplex()` construction
  - `is_freudenthal()` recognition function
  - `freudenthal_reflect()` structure-preserving reflection
  - Original bug regression tests
  - Adjacency and combinatorial properties

#### Integration and Edge Cases
- **`test_integration.jl`** - Integration tests and edge cases
  - Function composition workflows
  - Type consistency across all functions
  - Memory safety and mutation testing
  - Performance with large dimensions
  - Boundary conditions and degenerate cases
  - Comprehensive error handling
  - Numerical precision and stability

## Test Coverage

The test suite provides comprehensive coverage for all exported functions:

### Constructors and Properties
- All Simplex constructor variants
- Dimension and type functions
- Sorting and utility operations

### Geometric Operations
- General geometric reflection
- Mathematical property verification
- Error handling and edge cases

### Freudenthal Functions
- Initial simplex construction
- Freudenthal property recognition
- Structure-preserving reflection
- Integration with other functions

### Integration Tests
- Cross-function workflows
- Type safety and consistency
- Performance and scalability
- Edge cases and error conditions
