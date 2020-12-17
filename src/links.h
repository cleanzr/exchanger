#pragma once

#include "types.h"
#include <RcppArmadillo.h>
#include "records.h"
#include "entities.h"
#include "attribute_index.h"
#include "clust_params.h"

class Entities;
class Records;
class AbstractAttributeIndex;
class ClustParams;

/**
 * Represents the linkage stucture (links from entities to records)
 */
class Links {
private:
  std::vector<ent_id> forward_index_;

public:
  std::unordered_map<ent_id, std::unordered_set<rec_id>> inverse_index_;

public:
  /**
   * Setter for the linked entity of a record
   * @param rid a record identifier
   * @param eid an entity identifier
   */
  void set_link(rec_id rid, ent_id eid);

  /**
   * Getter for the linked entity of a record
   * @param rid a record identifier
   * @return the identifier of the linked entity
   */
  ent_id get_entity_id(rec_id rid) const;

  /**
   * Getter for the set of records linked to an entity
   * @param eid an entity identifier
   * @return a set of record identifiers
   */
  const std::unordered_set<rec_id>& get_record_ids(ent_id eid) const;

  /**
   * Convert to an R representation
   * @return a vector of entity assignments for the records
   */
  Rcpp::IntegerVector to_R() const;

  /**
   * Get the number of entities with linked records
   */
  int n_linked_ents() const;

  /**
   * Get the number of records linked to an entity
   * @param an entity identifier
   * @return the number of linked records
   */
  int n_records(ent_id eid) const;

  /**
   * Get the total number of records
   */
  int n_records() const;

  /**
   * Perform a Gibbs update
   * @param ents reference to entities container
   * @param recs reference to records container
   * @param clust_params reference to cluster parameters container
   * @param attr_indices a vector of attribute indices
   */
  void update(Entities &ents, const Records &recs, std::shared_ptr<ClustParams> clust_params, 
    const std::vector<std::shared_ptr<AbstractAttributeIndex> > &attr_indices);

private:
  /**
   * Perform a Gibbs update for the linked entity of a particular records
   * @param rid identifier of the record to update
   * @param ents reference to entities container
   * @param recs reference to records container
   * @param clust_params reference to cluster parameters container
   * @param attr_indices a vector of attribute indices
   */
  void update_link(rec_id rid, Entities &ents, const Records &recs, std::shared_ptr<ClustParams> clust_params,
    const std::vector<std::shared_ptr<AbstractAttributeIndex> > &attr_indices);

  /**
   * Add a link to the inverted index
   * @param rid a record identifier
   * @param eid an entity identifier
   */
  void add_inverted(rec_id rid, ent_id eid);

  /**
   * Remove a link from the inverted index
   * @param rid a record identifier
   * @param eid an entity identifier
   */
  void remove_inverted(rec_id rid, ent_id eid);
public:
  /**
   * Constructor from the R representation
   * @param forward_index a membership vector
   * @param clust_params RP S4 object representing the clustering parameters
   * @praam clust_prior RP S4 object representing the prior over the clustering parameters
   * @return a Links instance
   */
  Links(const Rcpp::IntegerVector &forward_index);
};
