# qSim Documentation

This directory contains the complete documentation for the qSim quantum computing simulation library.

## Document Overview

### 📋 [PRD.md](../PRD.md) - Product Requirements Document
**Audience**: Product managers, stakeholders, developers

**Purpose**: Defines what qSim will be, who it's for, and what problems it solves

**Contains**:
- Product vision and goals
- Target audience and use cases
- Feature requirements (Must-have, Should-have, Nice-to-have)
- Success metrics and KPIs
- Release criteria
- Out of scope items
- Future roadmap

**When to read**: Understanding project goals, planning features, making product decisions

---

### 🏗️ [ARCHITECTURE.md](./ARCHITECTURE.md) - System Architecture
**Audience**: Software architects, senior developers

**Purpose**: Explains the overall system design and architectural decisions

**Contains**:
- Layer architecture and module design
- Design principles and patterns
- Component interactions and data flow
- Concurrency model
- Platform abstraction strategy
- Extension points for future development
- Design decisions log

**When to read**: Understanding system structure, planning major changes, onboarding architects

---

### 🔧 [TECHNICAL_SPEC.md](./TECHNICAL_SPEC.md) - Technical Specification
**Audience**: Implementers, developers

**Purpose**: Detailed technical specifications for implementation

**Contains**:
- Data structure definitions and memory layout
- Algorithm implementations
- Gate application strategies
- BLAS integration details
- Performance optimization techniques
- Numerical precision considerations
- Platform-specific implementation details

**When to read**: Implementing features, debugging, optimizing performance

---

### 📚 [API.md](./API.md) - API Reference
**Audience**: Library users, application developers

**Purpose**: Complete reference for the public qSim API

**Contains**:
- Function signatures and parameters
- Usage examples
- Error handling patterns
- Code snippets for common tasks
- Complete API coverage

**When to read**: Using qSim in your application, learning the library, looking up function details

---

## Documentation Hierarchy

```
                    ┌─────────────┐
                    │   PRD.md    │  "What & Why"
                    │  (Product)  │
                    └──────┬──────┘
                           │
         ┌─────────────────┼─────────────────┐
         │                 │                 │
┌────────▼────────┐ ┌─────▼─────┐ ┌─────────▼────────┐
│ ARCHITECTURE.md │ │ TECHNICAL │ │     API.md       │
│   (Structure)   │ │  _SPEC.md │ │   (Interface)    │
│     "How"       │ │  (Details)│ │    "Usage"       │
└─────────────────┘ └───────────┘ └──────────────────┘
```

## Reading Guide

### For New Contributors
1. Start with [PRD.md](../PRD.md) to understand the project vision
2. Read [ARCHITECTURE.md](./ARCHITECTURE.md) to understand the system structure
3. Consult [API.md](./API.md) for public interfaces
4. Refer to [TECHNICAL_SPEC.md](./TECHNICAL_SPEC.md) when implementing features

### For Users of the Library
1. Read [API.md](./API.md) for usage instructions
2. Consult [PRD.md](../PRD.md) to understand supported features and roadmap
3. Check [TECHNICAL_SPEC.md](./TECHNICAL_SPEC.md) for performance characteristics

### For Maintainers
1. Keep [PRD.md](../PRD.md) updated with feature decisions
2. Update [ARCHITECTURE.md](./ARCHITECTURE.md) when making structural changes
3. Maintain [TECHNICAL_SPEC.md](./TECHNICAL_SPEC.md) as implementation evolves
4. Keep [API.md](./API.md) in sync with code changes

## Additional Documentation

### Source Code Documentation
- Inline comments in source files (`src/*.c`)
- Header documentation in public headers (`include/*.h`)

### Build Documentation
- Build instructions in root [README.md](../README.md)
- Platform-specific notes in `cmake/` files

### Test Documentation
- Test organization in `test/README.md` (to be created)
- Test utilities in `test/test.h`

## Document Status

| Document | Status | Last Updated | Completeness |
|----------|--------|--------------|--------------|
| PRD.md | Draft | 2025-12-22 | 100% (v1.0 scope) |
| ARCHITECTURE.md | Draft | 2025-12-22 | 80% (needs implementation decisions) |
| TECHNICAL_SPEC.md | Draft | 2025-12-22 | 60% (TBD sections) |
| API.md | Draft | 2025-12-22 | 40% (v1.0 APIs not implemented) |

## Contributing to Documentation

### Documentation Standards
- Use Markdown format
- Include table of contents for documents > 3 sections
- Add code examples for technical content
- Keep cross-references up to date
- Update "Last Updated" field when modifying

### Review Process
1. Documentation changes should accompany code changes
2. PRD changes require stakeholder approval
3. API documentation must match implementation
4. Architecture changes should be discussed before implementation

## Questions or Feedback?

- For product questions: Refer to [PRD.md](../PRD.md)
- For technical questions: Check [TECHNICAL_SPEC.md](./TECHNICAL_SPEC.md) and [API.md](./API.md)
- For architecture discussions: See [ARCHITECTURE.md](./ARCHITECTURE.md)
- For general questions: Create an issue in the project repository
