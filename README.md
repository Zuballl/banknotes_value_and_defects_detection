# Banknote Value and Defects Detection System

An automated vision-based system for identifying Polish banknote denominations and detecting mechanical and visual defects using MATLAB and computer vision techniques.

## Overview

This project implements a machine vision system that:
- Classifies banknotes by denomination (10, 20, 50, 100, 200 PLN)
- Performs geometric alignment relative to reference templates
- Detects and classifies surface defects (scratches, tears, creases, discoloration)
- Produces quality verdicts (accepted / requires inspection)

The system combines feature-based matching, image segmentation, and multi-channel defect analysis to provide robust denomination recognition and damage assessment.

## Requirements

- MATLAB R2020a or newer
- Computer Vision Toolbox
- Image Processing Toolbox

## Project Structure

```
├── setupTemplates.m          # Feature extraction from reference images
├── check_banknote.m          # Main verification script
├── templates/                # Reference banknote images
└── test/                     # Test images for validation
```

## Getting Started

### 1. Setup Template Features

Before running the verification system, extract features from reference images:

```matlab
run setupTemplates.m
```

This script processes the reference banknote images (one per denomination, three variants each) and saves extracted ORB features to `templateFeatures.mat`.

### 2. Run Banknote Verification

After setup, analyze test images:

```matlab
run check_banknote.m
```

The program will prompt for a test image path and perform:
- Denomination classification via feature matching
- Geometric registration to reference templates
- Defect detection and severity assessment

## Technical Approach

### Denomination Recognition
- **ORB (Oriented FAST and Rotated BRIEF)** – robust feature point detection and binary descriptors
- **Feature Matching** – Hamming distance-based descriptor comparison
- **Projective Transformation** – geometric alignment using `estimateGeometricTransform2D`

### Image Segmentation
- **Otsu Thresholding** – automatic binarization in HSV saturation channel
- **Morphological Operations** – closing, hole filling, and small region removal

### Mechanical Defect Detection
- **Mask Comparison** – logical difference between detected and reference masks
- **Connected Component Analysis** – identification of missing fragments
- **Edge-based Filtering** – artifact removal from registration errors

### Visual Anomaly Detection
- **LAB Color Space** – perceptually-uniform color difference analysis (ΔE metric)
- **Gaussian Smoothing** – color normalization before comparison
- **Adaptive Thresholding** – μ + 3σ criterion for color deviation detection

### Defect Classification
- **Region Properties** – area, eccentricity, edge proximity analysis
- **Two-Tier Thresholding** – minor defects (≤1.5%) vs. major defects (>1.5% or edge-located)

## Configuration Parameters

Key tunable parameters in `check_banknote.m`:

- `MAX_FEATURES` – maximum features to extract per image (default: 10,000)
- `MATCH_THRESHOLD` – feature matching similarity threshold (0.8)
- `MIN_DEFECT_AREA` – minimum detectable defect size in pixels (150)
- `MAJOR_AREA_PCT` – area threshold for major defect classification (1.5%)
- `EDGE_MARGIN` – pixel margin for edge region exclusion (12)
- `ACCEPTANCE.MINOR_MAX_COUNT` – maximum allowed minor defects (3)

## Output

The system displays:
- Annotated image with detected defects highlighted
- Individual defect reports (location, type, severity)
- Final verdict: **ACCEPTED** (no significant flaws) or **REQUIRES INSPECTION** (defects detected)


