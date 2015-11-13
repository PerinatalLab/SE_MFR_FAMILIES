SELECT grfar, grmor,
	GROUP_CONCAT(dads) as alldads,
	GROUP_CONCAT(kids) as allkids
FROM (
	SELECT p1.lpnr_barn as kids, p1.LopNrFar as dads,
		p2.LopNrFar as grfar, p2.LopNrMor as grmor
	FROM parents p1
	INNER JOIN parents p2
	ON p2.lpnr_barn=p1.LopNrFar
	WHERE p1.LopNrFar>0
) tt WHERE grfar>0 AND grmor>0
GROUP BY grfar, grmor
HAVING COUNT(DISTINCT dads)>1
INTO OUTFILE '/tmp/output_paternal.csv'
FIELDS TERMINATED BY ',' LINES TERMINATED BY '\n';

SELECT grfar, grmor,
        GROUP_CONCAT(moms) as allmoms,
        GROUP_CONCAT(kids) as allkids
FROM (
        SELECT p1.lpnr_barn as kids, p1.LopNrMor as moms,
                p2.LopNrFar as grfar, p2.LopNrMor as grmor
        FROM parents p1
        INNER JOIN parents p2
        ON p2.lpnr_barn=p1.LopNrMor
        WHERE p1.LopNrMor>0
) tt WHERE grfar>0 AND grmor>0
GROUP BY grfar, grmor
HAVING COUNT(DISTINCT moms)>1
INTO OUTFILE '/tmp/output_maternal.csv'
FIELDS TERMINATED BY ',' LINES TERMINATED BY '\n';


