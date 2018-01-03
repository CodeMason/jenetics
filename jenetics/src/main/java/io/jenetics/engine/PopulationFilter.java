package io.jenetics.engine;

import java.util.function.Predicate;
import java.util.function.Supplier;

import io.jenetics.Gene;
import io.jenetics.Phenotype;
import io.jenetics.util.Seq;

public interface PopulationFilter<
    G extends Gene<?, G>,
    C extends Comparable<? super C>> {
  FilterResult<G, C> filter(
      Predicate<? super Phenotype<G, C>> validator,
      Seq<Phenotype<G, C>> population,
      Supplier<Phenotype<G, C>> phenotypeSupplier,
      long maximalPhenotypeAge,
      long generation
  );
}
