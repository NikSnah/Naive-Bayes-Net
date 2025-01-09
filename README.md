# Naive-Bayes-Net
A UOFT csc384 project that implements a Naive Bayes Network, a probabilistic graphical model based on the Bayesian Theorem

---

## Features

- **Bayesian Framework**:
  - Uses **prior probabilities** and **likelihood** to compute posterior probabilities.
- **Naive Bayes Assumption**:
  - Assumes conditional independence between features for computational efficiency.
- **Custom Implementation**:
  - Modular design for flexibility and reusability.
- **Autograder**:
  - Automated testing of the implementation for correctness and robustness.

---

## File Overview

### **1. `bnetbase.py`**
- Provides the core implementation of a **Bayes Network**.
- Handles:
  - Node creation and management.
  - Conditional probability tables (CPTs).
  - Inference calculations using Bayesian principles.

### **2. `naive_bayes_solution.py`**
- Implements the **Naive Bayes algorithm** for classification.
- Includes methods for:
  - Training the model on labeled datasets.
  - Predicting class probabilities for unseen data.
  - Evaluating the model's accuracy.
- Demonstrates practical applications of Naive Bayes in real-world scenarios.

### **3. `autograder.py`**
- An automated testing script for validating the implementation.
- Evaluates:
  - Correctness of probabilistic computations.
  - Accuracy of predictions against known outputs.
  - Efficiency of the code.

---

## How It Works

1. **Model Initialization**:
   - Define features and class labels.
   - Provide prior probabilities and likelihood functions.

2. **Training**:
   - Use labeled data to compute probabilities and build the model.

3. **Prediction**:
   - Input new data to compute posterior probabilities and classify.

4. **Evaluation**:
   - Test the model's accuracy using the autograder or custom test cases.

---

## Installation

1. **Clone the Repository**:
   ```bash
   git clone https://github.com/username/naive-bayes-network.git
   cd naive-bayes-network
