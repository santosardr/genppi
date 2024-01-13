CREATE OR REPLACE FUNCTION count_intersect()
RETURNS TABLE(genome text, query_nounique bigint, query_unique bigint, string_intersect_unique bigint, unique_intersect_ratio decimal, string_proportion decimal) AS $$
BEGIN
    RETURN QUERY
WITH stringset AS (
  SELECT DISTINCT 
    LEAST(id1, id2) AS id1, 
    GREATEST(id1, id2) AS id2
  FROM related 
  WHERE genome1 = 'string'
)    
SELECT 
  final.genome, 
  final.query_nounique, 
  final.query_unique,  
  final.string_intersect_unique,
  round(final.query_unique::DECIMAL / final.string_intersect_unique,2) AS unique_intersect_ratio,
  round(final.string_intersect_unique::DECIMAL / final.string_interactions, 2) AS string_proportion
FROM (
  SELECT 
    query.genome1::text AS genome, 
    count(1) AS query_nounique, 
    (
      SELECT count(1)
      FROM (
	      SELECT LEAST(id1, id2), GREATEST(id1, id2)
    	  FROM(
		        SELECT 
				  related.id1, 
				  related.id2
    		    FROM related, 
    		    stringset 
    		    WHERE genome1 = query.genome1
    		      AND    LEAST(related.id1, related.id2) = stringset.id1 
    		      AND GREATEST(related.id1, related.id2) = stringset.id2
    		  ) AS into_string
    	  ) as level2 
    ) AS string_intersect_unique,
    (SELECT count(*) FROM (
    	SELECT DISTINCT LEAST(id1, id2), GREATEST(id1, id2)
    	FROM related
    	WHERE genome1 = query.genome1
    ) AS distinct_pairs) AS query_unique,
    (select count(1) from stringset) as string_interactions
  FROM related AS query 
  GROUP BY query.genome1
  ORDER BY string_intersect_unique
) AS final 
ORDER BY unique_intersect_ratio;
END;
$$ LANGUAGE plpgsql;
select * from count_intersect();
