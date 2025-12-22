# State API Documentation

**Module:** State Representation
**Header:** `include/state.h`
**Implementation:** `src/state.c`

---

## Data Types

### `state_type_t`

```c
typedef enum state_type {
    PURE,   // Pure quantum state (state vector)
    MIXED   // Mixed quantum state (density matrix, lower triangle)
} state_type_t;
```

**Description:** Enumeration distinguishing pure and mixed quantum states.

---

### `state_t`

```c
typedef struct state {
    state_type_t type;  // State type (PURE or MIXED)
    cplx_t *data;       // Complex array storing state representation
    qubit_t qubits;     // Number of qubits
} state_t;
```

**Description:** Structure representing a multi-qubit quantum state.

**Fields:**
- `type`: Determines representation format (state vector vs density matrix)
- `data`: Pointer to complex array containing state information
  - **PURE:** Array of length `2^qubits` (state vector)
  - **MIXED:** Array of length `2^qubits * (2^qubits + 1) / 2` (lower triangle of density matrix)
- `qubits`: Number of qubits in the system

**Memory Ownership:**
- The `state_t` struct owns the memory pointed to by `data`
- User must call `state_free()` to release memory
- Copying a `state_t` requires `state_cp()` for deep copy; shallow copy causes double-free

---

## Functions

### `state_free()`

```c
void state_free(state_t *state);
```

**Description:** Deallocates memory for the state representation and resets metadata.

**Parameters:**
- `state` (in/out): Pointer to state to be freed

**Behavior:**
- Frees `state->data` and sets it to `NULL`
- Resets `state->qubits` to 0
- Safe to call multiple times on the same state
- Does NOT free the `state_t` struct itself (only the `data` array)

**Example:**
```c
state_t state = {.type = PURE, .data = NULL, .qubits = 0};
state_init(&state, 3, NULL);
// ... use state ...
state_free(&state);  // Deallocates state.data
```

---

### `state_len()`

```c
dim_t state_len(const state_t *state);
```

**Description:** Computes the number of complex elements in the state representation.

**Parameters:**
- `state` (in): Pointer to state (only `type` and `qubits` fields are accessed)

**Returns:**
- **PURE state:** `2^qubits` (dimension of Hilbert space)
- **MIXED state:** `2^qubits * (2^qubits + 1) / 2` (lower triangle of density matrix)
- Returns 0 if `state->qubits` is 0

**Example:**
```c
state_t state = {.type = PURE, .qubits = 3};
dim_t len = state_len(&state);  // Returns 8

state.type = MIXED;
len = state_len(&state);  // Returns 36
```

**Note:** Does NOT validate that `state->data` is allocated; only computes length based on metadata.

---

### `state_init()`

```c
void state_init(state_t *state, qubit_t qubits, cplx_t **data);
```

**Description:** Initializes a quantum state with the specified number of qubits.

**Parameters:**
- `state` (in/out): Pointer to state structure to initialize
  - `state->type` must be set before calling (determines array length)
- `qubits` (in): Number of qubits for the system
- `data` (in/out): Pointer to pointer to complex array (double pointer for ownership transfer)
  - If `*data != NULL`: Takes ownership of the array (sets `*data = NULL` after transfer)
  - If `*data == NULL`: Allocates and zero-initializes array

**Behavior:**
- If `state->data` is already allocated, frees old memory first (safe for reinitialization)
- Sets `state->qubits = qubits`
- Either takes ownership of provided array or allocates new zero-initialized array
- Prints error to `stderr` and sets `state->qubits = 0` if allocation fails

**Memory Ownership:**
- If `data != NULL` and `*data != NULL`, the function **steals** ownership of `*data`
- After the call, `*data == NULL` (prevents double-free by caller)
- Caller must NOT free the array after passing it to `state_init()`

**Example 1: Zero initialization**
```c
state_t state = {.type = PURE, .data = NULL, .qubits = 0};
state_init(&state, 2, NULL);
// state.data now points to zero-initialized array of length 4
```

**Example 2: Ownership transfer**
```c
state_t state = {.type = PURE, .data = NULL, .qubits = 0};

cplx_t *my_data = malloc(4 * sizeof(cplx_t));
my_data[0] = 1.0; my_data[1] = 0.0; my_data[2] = 0.0; my_data[3] = 0.0;

state_init(&state, 2, &my_data);
// my_data is now NULL (ownership transferred)
// state.data points to the array (don't free my_data!)
```

**Example 3: Reinitialization**
```c
state_t state = {.type = PURE, .data = NULL, .qubits = 0};
state_init(&state, 2, NULL);  // Allocates 4 elements

state_init(&state, 3, NULL);  // Frees old array, allocates 8 elements
```

---

### `state_plus()`

```c
void state_plus(state_t *state, qubit_t qubits);
```

**Description:** Initializes a quantum state to the uniform superposition of all computational basis states (plus state).

**Parameters:**
- `state` (in/out): Pointer to state structure
  - `state->type` must be set before calling
- `qubits` (in): Number of qubits

**Behavior:**
- **PURE state:** Sets all amplitudes to `1/sqrt(2^qubits)` (normalized uniform superposition)
  - Represents: `|+⟩^⊗n = (|0⟩ + |1⟩)^⊗n / sqrt(2^n)`
- **MIXED state:** Sets all density matrix entries to `1/2^qubits`
  - Represents: `ρ = I / 2^n` (maximally mixed state)
- Calls `state_init()` internally (handles memory allocation)
- If allocation fails, frees state and prints error to `stderr`

**Example:**
```c
// Pure state: uniform superposition
state_t psi = {.type = PURE, .data = NULL, .qubits = 0};
state_plus(&psi, 2);
// psi.data = [0.5, 0.5, 0.5, 0.5] (all real, normalized)

// Mixed state: maximally mixed
state_t rho = {.type = MIXED, .data = NULL, .qubits = 0};
state_plus(&rho, 2);
// rho represents I/4 (all 10 lower-triangle entries = 0.25)
```

**Note:** This is a convenience function. For custom initial states, use `state_init()` directly.

---

### `state_cp()`

```c
state_t state_cp(const state_t *state);
```

**Description:** Creates a deep copy of a quantum state.

**Parameters:**
- `state` (in): Pointer to state to copy

**Returns:**
- Deep copy of the state (new allocation for `data` array)
- If allocation fails, returns state with `data = NULL` and `qubits = 0`, prints error to `stderr`

**Behavior:**
- Copies all metadata (`type`, `qubits`)
- Allocates new array and copies all elements using `cblas_zcopy()`
- Original and copy are completely independent (modifying one does not affect the other)

**Memory Management:**
- Returned state must be freed with `state_free()` when done
- Does NOT free the original state

**Example:**
```c
state_t original = {.type = PURE, .data = NULL, .qubits = 0};
state_plus(&original, 3);

state_t copy = state_cp(&original);
// copy.data != original.data (different allocations)

original.data[0] = 999.0;  // Does NOT affect copy

state_free(&original);
state_free(&copy);  // Must free both independently
```

**Error Handling:**
- If allocation fails, prints `"state_cp(): out.data allocation failed"` to `stderr`
- Returns a state with `data = NULL`, `qubits = 0`
- Caller should check `copy.data != NULL` before use

---

## Usage Patterns

### Pattern 1: Create and Initialize

```c
state_t state = {.type = PURE, .data = NULL, .qubits = 0};
state_init(&state, num_qubits, NULL);

// ... use state ...

state_free(&state);
```

### Pattern 2: Initialize with Custom Data

```c
state_t state = {.type = PURE, .data = NULL, .qubits = 0};

cplx_t *data = malloc(4 * sizeof(cplx_t));
// ... populate data ...

state_init(&state, 2, &data);  // Transfers ownership
// data is now NULL, don't free it!

state_free(&state);  // Frees the transferred array
```

### Pattern 3: Uniform Superposition

```c
state_t state = {.type = PURE, .data = NULL, .qubits = 0};
state_plus(&state, 3);

// ... use state ...

state_free(&state);
```

### Pattern 4: Deep Copy

```c
state_t original = {.type = MIXED, .data = NULL, .qubits = 0};
state_init(&original, 2, NULL);

state_t copy = state_cp(&original);

// Modify copy independently
copy.data[0] = 1.0;

state_free(&original);
state_free(&copy);
```

---

## Common Pitfalls

### ❌ Shallow Copy (Double-Free)
```c
state_t state1 = {.type = PURE, .data = NULL, .qubits = 0};
state_init(&state1, 2, NULL);

state_t state2 = state1;  // WRONG: Shallow copy!

state_free(&state1);
state_free(&state2);  // Double-free! Undefined behavior
```

**✅ Correct: Deep Copy**
```c
state_t state2 = state_cp(&state1);
state_free(&state1);
state_free(&state2);  // Safe
```

---

### ❌ Using Data After Ownership Transfer
```c
cplx_t *data = malloc(4 * sizeof(cplx_t));
state_init(&state, 2, &data);

data[0] = 1.0;  // WRONG: data is NULL now!
```

**✅ Correct: Access via state.data**
```c
cplx_t *data = malloc(4 * sizeof(cplx_t));
state_init(&state, 2, &data);

state.data[0] = 1.0;  // Correct
```

---

### ❌ Forgetting to Set type
```c
state_t state;  // Uninitialized!
state_init(&state, 2, NULL);  // state.type is garbage
```

**✅ Correct: Initialize struct**
```c
state_t state = {.type = PURE, .data = NULL, .qubits = 0};
state_init(&state, 2, NULL);
```

---

## Thread Safety

**None of these functions are thread-safe.**

- Multiple threads can read the same `state_t` concurrently (read-only operations)
- Modifying a `state_t` from multiple threads requires external synchronization
- Allocating separate states in different threads is safe

---

## Error Handling

All functions use **error-by-convention**:
- Allocation failures print to `stderr`
- Invalid states have `data = NULL` and `qubits = 0`
- No exceptions or error codes returned

**Checking for errors:**
```c
state_t state = {.type = PURE, .data = NULL, .qubits = 0};
state_init(&state, 10, NULL);

if (state.data == NULL) {
    // Allocation failed
    fprintf(stderr, "Failed to initialize state\n");
    return -1;
}
```
