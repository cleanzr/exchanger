# exchanger 0.3.0
* Fixes error in Gibbs update for entity attributes and links
* Uses improved scaling for "temperature" of distortion softmax distribution
* Adds support for point mass priors on the distortion probabilities
* Adds support for a gamma prior on the concentration parameter associated with 
the distortion distribution
* Fixes bug in update for Pitman-Yor parameters

# exchanger 0.2.0
* First release after major refactoring
* Depend on `comparator` and `clevr` for distance functions and evaluation 
  functions respectively
* Use underscore_case instead of camelCase for function/variable names
* Implement show methods for user-facing objects
* Improve quality of documentation