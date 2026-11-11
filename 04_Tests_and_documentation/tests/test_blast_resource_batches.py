"""Bound BLAST's internal query grouping without changing search parameters.

The installed BLAST reader consumes complete queries until their total bases
meet or exceed BATCH_SIZE. Queries are ordered arm, full-primary per candidate.
These tests independently reconstruct that threshold loop: a resource grouping
must never combine queries belonging to different candidates.
"""
import copy
import random
import sys
import unittest
from pathlib import Path

sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
import design_engine as engine


def candidates(lengths,*,prefix='CCGATTCGTTAGTT',suffix='TTCCGATTCGTTAG'):
    result=[]
    for index,length in enumerate(lengths):
        arm=('ACGT'*((length+3)//4))[:length]
        result.append(dict(id=f'candidate_{index}',arm=arm,full=prefix+arm+suffix))
    return result


def threshold_groups(items,threshold):
    """Independent reconstruction of NCBI's total-bases-until-threshold loop."""
    pending=[];bases=0;groups=[]
    for candidate in items:
        for kind in ('arm','full'):
            pending.append((candidate['id'],kind));bases+=len(candidate[kind])
            if bases>=threshold:
                groups.append(pending);pending=[];bases=0
    if pending:groups.append(pending)
    return groups


class BlastResourceBatchTests(unittest.TestCase):
    def assert_no_cross_candidate_group(self,items,limit):
        groups=threshold_groups(items,limit)
        expected=[(c['id'],kind) for c in items for kind in ('arm','full')]
        self.assertEqual([query for group in groups for query in group],expected)
        self.assertTrue(all(len({candidate for candidate,_ in group})==1 for group in groups))
        return groups

    def test_current_errfi1_36_through_90_arms_keep_each_pair_together_at_91(self):
        items=candidates(range(36,91))
        self.assertGreaterEqual(min(len(c['arm'])+len(c['full']) for c in items),100)
        limit=engine._blast_internal_batch_letters(items)
        self.assertEqual(limit,91)
        groups=self.assert_no_cross_candidate_group(items,limit)
        self.assertEqual(len(groups),len(items))
        self.assertTrue(all([kind for _,kind in group]==['arm','full'] for group in groups))

    def test_new_docker_lengths_and_spacers_do_not_change_pair_ownership(self):
        dockers=['tattcgtTCCTCCTCCTCCTCCTCCT','tattcgtCAACAACAACAACAACAAC',
                 'tattcgtACCACCACCACCACCACCA','tattcgtCTCTCTCTCTCTCTCTCTC',
                 'tattcgtAACAACAACAACAACAACAA']
        items=[]
        for index,(length,docker) in enumerate(zip([36,42,57,76,90],dockers)):
            item=candidates([length],suffix='TT'+docker)[0]
            item['id']=f'docker_candidate_{index}';items.append(item)
        before=copy.deepcopy(items)
        limit=engine._blast_internal_batch_letters(items)
        self.assertEqual(limit,91)
        self.assertEqual(len(self.assert_no_cross_candidate_group(items,limit)),len(items))
        self.assertEqual(items,before,'Resource grouping must not edit any sequence or candidate')

    def test_threshold_equality_finishes_the_pair_before_reading_next_candidate(self):
        # Largest arm is 99; smallest arm/full pair has exactly 100 bases.
        items=candidates([36,99,40])
        self.assertEqual(len(items[0]['arm'])+len(items[0]['full']),100)
        self.assertEqual(engine._blast_internal_batch_letters(items),100)
        self.assertEqual(len(self.assert_no_cross_candidate_group(items,100)),3)

    def test_wide_arm_length_mixture_falls_back_to_one_complete_query_at_a_time(self):
        items=candidates([36,180,42])
        limit=engine._blast_internal_batch_letters(items)
        self.assertEqual(limit,1)
        groups=self.assert_no_cross_candidate_group(items,limit)
        self.assertEqual(len(groups),2*len(items))
        self.assertTrue(all(len(group)==1 for group in groups))

    def test_order_and_trailing_partial_outer_groups_keep_all_queries(self):
        for lengths in ([90,36,42,60],[36],[90],[90,36],[42,90,36]):
            with self.subTest(lengths=lengths):
                items=candidates(lengths)
                limit=engine._blast_internal_batch_letters(items)
                self.assertEqual(limit,max(lengths)+1)
                self.assertEqual(len(self.assert_no_cross_candidate_group(items,limit)),len(items))

    def test_brute_varied_positive_lengths_never_share_internal_candidate_group(self):
        rng=random.Random(20261103)
        for trial in range(1000):
            items=[]
            for index in range(rng.randrange(1,9)):
                length=rng.randrange(1,251)
                item=candidates([length],prefix='A'*rng.randrange(0,51),suffix='T'*rng.randrange(0,51))[0]
                item['id']=f'{trial}_{index}';items.append(item)
            limit=engine._blast_internal_batch_letters(items)
            pair_limit=max(len(c['arm']) for c in items)+1
            expected=pair_limit if pair_limit<=min(len(c['arm'])+len(c['full']) for c in items) else 1
            self.assertEqual(limit,expected)
            self.assert_no_cross_candidate_group(items,limit)


if __name__=='__main__':unittest.main()
