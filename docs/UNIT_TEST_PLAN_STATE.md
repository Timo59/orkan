# Unit Test Plan: state.c

## Test Structure

Each function should be tested with:
- **PURE** states (separate test)
- **MIXED** states (separate test)

---

## 1. `state_len()` Tests

### Test: `test_state_len_pure()`
**Purpose:** Verify length calculation for pure states

| Qubits | Expected Length | Formula |
|--------|----------------|---------|
| 1      | 2              | 2^1 |
| 2      | 4              | 2^2 |
| 3      | 8              | 2^3 |
| 4      | 16             | 2^4 |

**Implementation:**
```c
void test_state_len_pure(void) {
    state_t state;
    state.type = PURE;

    state.qubits = 1;
    TEST_ASSERT_EQUAL_INT64(2, state_len(&state));

    state.qubits = 2;
    TEST_ASSERT_EQUAL_INT64(4, state_len(&state));

    state.qubits = 3;
    TEST_ASSERT_EQUAL_INT64(8, state_len(&state));

    state.qubits = 4;
    TEST_ASSERT_EQUAL_INT64(16, state_len(&state));
}
```

### Test: `test_state_len_mixed()`
**Purpose:** Verify length calculation for mixed states (lower triangle)

| Qubits | Dimension | Expected Length | Formula |
|--------|-----------|----------------|---------|
| 1      | 2         | 3              | 2×3/2 |
| 2      | 4         | 10             | 4×5/2 |
| 3      | 8         | 36             | 8×9/2 |

**Implementation:**
```c
void test_state_len_mixed(void) {
    state_t state;
    state.type = MIXED;

    state.qubits = 1;
    TEST_ASSERT_EQUAL_INT64(3, state_len(&state));

    state.qubits = 2;
    TEST_ASSERT_EQUAL_INT64(10, state_len(&state));

    state.qubits = 3;
    TEST_ASSERT_EQUAL_INT64(36, state_len(&state));
}
```

---

## 2. `state_init()` Tests

### Test: `test_state_init_pure_null_data()`
**Purpose:** Initialize pure state with NULL data (zero initialization)

**Steps:**
1. Initialize state with 2 qubits, NULL data
2. Verify `state.qubits == 2`
3. Verify `state.data != NULL`
4. Verify all elements are zero
5. Call `state_free()` and verify no leaks

### Test: `test_state_init_pure_with_data()`
**Purpose:** Initialize pure state with provided data (ownership transfer)

**Steps:**
1. Allocate array with known values
2. Call `state_init()` with pointer to array
3. Verify array pointer is now NULL (ownership transferred)
4. Verify state.data points to the original allocation
5. Verify values are preserved
6. Call `state_free()`

### Test: `test_state_init_mixed_null_data()`
**Purpose:** Same as pure, but for MIXED state

### Test: `test_state_init_mixed_with_data()`
**Purpose:** Same as pure, but for MIXED state

### Test: `test_state_init_reinitialization()`
**Purpose:** Verify that reinitializing a state frees old memory

**Steps:**
1. Initialize state with 2 qubits
2. Store pointer to data
3. Re-initialize same state with 3 qubits
4. Verify new allocation occurred
5. Verify old pointer is invalid (hard to test without ASAN)
6. Verify state.qubits == 3

**Code example:**
```c
void test_state_init_pure_null_data(void) {
    state_t state = {.type = PURE, .data = NULL, .qubits = 0};

    state_init(&state, 2, NULL);

    TEST_ASSERT_EQUAL_UINT8(2, state.qubits);
    TEST_ASSERT_NOT_NULL(state.data);

    // Verify zero initialization
    dim_t len = state_len(&state);
    for (dim_t i = 0; i < len; i++) {
        TEST_ASSERT_COMPLEX_WITHIN(0.0 + 0.0*I, state.data[i], PRECISION);
    }

    state_free(&state);
}
```

```c
void test_state_init_pure_with_data(void) {
    state_t state = {.type = PURE, .data = NULL, .qubits = 0};

    // Allocate and initialize data
    cplx_t *data = malloc(4 * sizeof(cplx_t));
    TEST_ASSERT_NOT_NULL(data);
    data[0] = 1.0 + 0.0*I;
    data[1] = 0.0 + 1.0*I;
    data[2] = -1.0 + 0.0*I;
    data[3] = 0.0 - 1.0*I;

    cplx_t *original_ptr = data;

    state_init(&state, 2, &data);

    // Verify ownership transfer
    TEST_ASSERT_NULL(data);
    TEST_ASSERT_EQUAL_PTR(original_ptr, state.data);

    // Verify values preserved
    TEST_ASSERT_COMPLEX_WITHIN(1.0 + 0.0*I, state.data[0], PRECISION);
    TEST_ASSERT_COMPLEX_WITHIN(0.0 + 1.0*I, state.data[1], PRECISION);
    TEST_ASSERT_COMPLEX_WITHIN(-1.0 + 0.0*I, state.data[2], PRECISION);
    TEST_ASSERT_COMPLEX_WITHIN(0.0 - 1.0*I, state.data[3], PRECISION);

    state_free(&state);
}
```

---

## 3. `state_free()` Tests

### Test: `test_state_free_pure()`
**Purpose:** Verify memory deallocation and reset

**Steps:**
1. Initialize state with 3 qubits
2. Call `state_free()`
3. Verify `state.data == NULL`
4. Verify `state.qubits == 0`
5. Verify can be called again safely (double-free safety)

### Test: `test_state_free_mixed()`
**Purpose:** Same for MIXED state

**Code example:**
```c
void test_state_free_pure(void) {
    state_t state = {.type = PURE, .data = NULL, .qubits = 0};

    state_init(&state, 3, NULL);
    TEST_ASSERT_NOT_NULL(state.data);
    TEST_ASSERT_EQUAL_UINT8(3, state.qubits);

    state_free(&state);

    TEST_ASSERT_NULL(state.data);
    TEST_ASSERT_EQUAL_UINT8(0, state.qubits);

    // Safe to call again
    state_free(&state);
    TEST_ASSERT_NULL(state.data);
}
```

---

## 4. `state_plus()` Tests

### Test: `test_state_plus_pure_normalization()`
**Purpose:** Verify uniform superposition is normalized

**Steps:**
1. Create plus state with 3 qubits
2. Verify all amplitudes equal `1/sqrt(8)`
3. Verify normalization: `sum |amplitude|^2 == 1`

### Test: `test_state_plus_mixed_normalization()`
**Purpose:** Verify mixed plus state has correct entries

**Steps:**
1. Create mixed plus state with 2 qubits
2. Verify diagonal entries equal `1/4`
3. Verify off-diagonal entries in lower triangle equal `1/4`
4. Verify trace == 1 (sum of diagonal)

**Code example:**
```c
void test_state_plus_pure_normalization(void) {
    state_t state = {.type = PURE, .data = NULL, .qubits = 0};

    state_plus(&state, 3);

    dim_t len = state_len(&state);
    TEST_ASSERT_EQUAL_INT64(8, len);

    double expected_amplitude = 1.0 / sqrt(8.0);

    // Check all amplitudes
    for (dim_t i = 0; i < len; i++) {
        double real_part = creal(state.data[i]);
        double imag_part = cimag(state.data[i]);
        TEST_ASSERT_DOUBLE_WITHIN(PRECISION, expected_amplitude, real_part);
        TEST_ASSERT_DOUBLE_WITHIN(PRECISION, 0.0, imag_part);
    }

    // Verify normalization
    double norm_sq = 0.0;
    for (dim_t i = 0; i < len; i++) {
        double mag = cabs(state.data[i]);
        norm_sq += mag * mag;
    }
    TEST_ASSERT_DOUBLE_WITHIN(PRECISION, 1.0, norm_sq);

    state_free(&state);
}
```

```c
void test_state_plus_mixed_trace(void) {
    state_t state = {.type = MIXED, .data = NULL, .qubits = 0};

    state_plus(&state, 2);

    dim_t dim = POW2(2, dim_t);  // Dimension = 4
    double expected_entry = 1.0 / (double)dim;  // 1/4

    // For mixed state, all entries should be 1/4
    dim_t len = state_len(&state);
    for (dim_t i = 0; i < len; i++) {
        double real_part = creal(state.data[i]);
        double imag_part = cimag(state.data[i]);
        TEST_ASSERT_DOUBLE_WITHIN(PRECISION, expected_entry, real_part);
        TEST_ASSERT_DOUBLE_WITHIN(PRECISION, 0.0, imag_part);
    }

    // Verify trace = 1 (sum of diagonal elements in packed storage)
    // Diagonal indices in packed storage: 0, 2, 5, 9 for 4x4 matrix
    // General formula: i*(i+1)/2 + i = i*(i+3)/2
    double trace = 0.0;
    for (dim_t i = 0; i < dim; i++) {
        dim_t idx = i * (i + 3) / 2;
        trace += creal(state.data[idx]);
    }
    TEST_ASSERT_DOUBLE_WITHIN(PRECISION, 1.0, trace);

    state_free(&state);
}
```

---

## 5. `state_cp()` Tests

### Test: `test_state_cp_pure_deep_copy()`
**Purpose:** Verify deep copy creates independent state

**Steps:**
1. Create original state with specific values
2. Create copy using `state_cp()`
3. Verify `copy.data != original.data` (different pointers)
4. Verify all values are identical
5. Modify original
6. Verify copy is unchanged
7. Free both states

### Test: `test_state_cp_mixed_deep_copy()`
**Purpose:** Same for MIXED state

**Code example:**
```c
void test_state_cp_pure_deep_copy(void) {
    state_t original = {.type = PURE, .data = NULL, .qubits = 0};

    state_plus(&original, 2);

    state_t copy = state_cp(&original);

    // Verify metadata copied
    TEST_ASSERT_EQUAL_UINT8(original.type, copy.type);
    TEST_ASSERT_EQUAL_UINT8(original.qubits, copy.qubits);

    // Verify different allocations
    TEST_ASSERT_NOT_EQUAL_PTR(original.data, copy.data);

    // Verify values identical
    dim_t len = state_len(&original);
    for (dim_t i = 0; i < len; i++) {
        TEST_ASSERT_COMPLEX_WITHIN(original.data[i], copy.data[i], PRECISION);
    }

    // Modify original
    original.data[0] = 99.0 + 99.0*I;

    // Verify copy unchanged
    double expected = 1.0 / sqrt(4.0);
    TEST_ASSERT_DOUBLE_WITHIN(PRECISION, expected, creal(copy.data[0]));

    state_free(&original);
    state_free(&copy);
}
```

---

## 6. Edge Case Tests

### Test: `test_state_single_qubit()`
**Purpose:** Verify 1-qubit states work correctly

### Test: `test_state_many_qubits()`
**Purpose:** Test with larger qubit counts (e.g., 10 qubits)
- Mainly to verify no integer overflow in `state_len()`
- Verify allocations succeed (or fail gracefully if too large)

---

## Summary: Test Coverage

| Function | # Tests | Coverage |
|----------|---------|----------|
| `state_len()` | 2 | PURE, MIXED with multiple qubit counts |
| `state_init()` | 5 | NULL data, with data, reinitialization for both types |
| `state_free()` | 2 | PURE, MIXED, double-free safety |
| `state_plus()` | 2 | PURE normalization, MIXED trace |
| `state_cp()` | 2 | Deep copy independence for both types |
| Edge cases | 2+ | Single qubit, many qubits, allocation limits |

**Total:** ~15-18 unit tests

---

## Notes on Implementation

1. **Use existing test infrastructure** from `test/include/test.h`:
   - `TEST_ASSERT_COMPLEX_WITHIN()` for complex comparisons
   - `PRECISION` macro for tolerance

2. **Memory leak detection**: Run tests under Valgrind or ASAN to verify no leaks

3. **Test organization**: Group by function in `test_state.c`:
   ```c
   void test_state_len_pure(void);
   void test_state_len_mixed(void);

   void test_state_init_pure_null_data(void);
   void test_state_init_pure_with_data(void);
   // ... etc
   ```

4. **Future additions**:
   - Once gate operations are implemented, add integration tests
   - Test state behavior after gate applications
