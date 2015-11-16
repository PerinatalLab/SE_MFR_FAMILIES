SELECT grfar, grmor,
	GROUP_CONCAT(moms) as allmoms,
	GROUP_CONCAT(kids) as allkids,
	GROUP_CONCAT(dads) as alldads,
	GROUP_CONCAT(birthyear) as birthyears,
	GROUP_CONCAT(dadage) as dadages,
	GROUP_CONCAT(momage) as momages
FROM (
	SELECT p1.lpnr_barn as kids, p1.LopNrMor as moms, p1.LopNrFar as dads,
		p2.LopNrFar as grfar, p2.LopNrMor as grmor,
		p1.FoddArBarn as birthyear,
		p1.FoddArBarn-p1.FoddArFar as dadage,
		p1.FoddArBarn-p1.FoddArMor as momage
	FROM parents p1
	INNER JOIN parents p2
	ON p2.lpnr_barn=p1.LopNrMor
	WHERE p1.LopNrMor>0
) tt WHERE grfar>0 AND grmor>0
GROUP BY grfar, grmor
HAVING COUNT(DISTINCT moms)>1
INTO OUTFILE '/tmp/output.csv' FIELDS TERMINATED BY '\t' LINES TERMINATED BY '\n';

/* finding cousins whose fathers are fullsibs (brothers) */
SELECT grfar, grmor,
	GROUP_CONCAT(dads) as alldads,
	GROUP_CONCAT(kids) as allkids,
	GROUP_CONCAT(moms) as allmoms,
	GROUP_CONCAT(birthyear) as birthyears,
	GROUP_CONCAT(dadage) as dadages,
	GROUP_CONCAT(momage) as momages
FROM (
	SELECT p1.lpnr_barn as kids, p1.LopNrMor as moms, p1.LopNrFar as dads,
		p2.LopNrFar as grfar, p2.LopNrMor as grmor,
		p1.FoddArBarn as birthyear,
		p1.FoddArBarn-p1.FoddArFar as dadage,
		p1.FoddArBarn-p1.FoddArMor as momage
	FROM parents p1
	INNER JOIN parents p2
	ON p2.lpnr_barn=p1.LopNrFar
	WHERE p1.LopNrFar>0
) tt WHERE grfar>0 AND grmor>0
GROUP BY grfar, grmor
HAVING COUNT(DISTINCT dads)>1
INTO OUTFILE '/tmp/output.csv' FIELDS TERMINATED BY '\t' LINES TERMINATED BY '\n';


