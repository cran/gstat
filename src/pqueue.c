/* 
 * priority queue, Jan-Apr 1998
 *
 * I did a net search one day, and found so many implementations of
 * priority queues that I decided to write my own. For the fun of it.
 *
 * This one has no limits and needs a generic (qsort-like) comparison
 * function for element comparison at initialisation. To enqueue an
 * unordered array of elements, elements are sorted with qsort(), before
 * they are merged into the ordered queue. Only pointers to the next
 * element are stored (so it's basically an ordered single linked list).
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "defs.h"
#include "userio.h"
#include "data.h"
#include "utils.h"
#include "nsearch.h"
#include "pqueue.h"

static void enlarge_queue(QUEUE *q);

static void enlarge_queue(QUEUE *q) {
/*
 * don't use realloc() on `empty': all pointers to the
 * first element would have to be moved along with itself.
 */
	int i;
	Q_ELEMENT *block;

	block = (Q_ELEMENT *) emalloc(Q_BUFFER_SIZE * sizeof(Q_ELEMENT));
	for (i = 0; i < Q_BUFFER_SIZE - 1; i++)
		block[i].next = &(block[i+1]);
	block[Q_BUFFER_SIZE - 1].next = NULL;

	if (q->empty == NULL)
		q->empty = block;
	else
		q->empty->next = block;
	q->max_length += Q_BUFFER_SIZE;

	/* register block to free later on... */
	q->blocks += 1;
	q->block = (Q_ELEMENT **) erealloc(q->block, q->blocks * sizeof(Q_ELEMENT *));
	q->block[q->blocks - 1] = block;
}

QUEUE *init_queue(QUEUE *q, int (CDECL *cmp)(const QUEUE_NODE *a, 
		const QUEUE_NODE *b)) {
	int i, j;

	if (q == NULL) {
		q = (QUEUE *) emalloc(sizeof(QUEUE));
		q->max_length = q->blocks = 0;
		q->empty = NULL;
		q->block = NULL;
		q->cmp = cmp;
		enlarge_queue(q);
	} else {
		q->empty = q->block[0];
		for (i = 0; i < q->blocks; i++) {
			for (j = 0; j < Q_BUFFER_SIZE - 1; j++) /* connect elements: */
				q->block[i][j].next = &(q->block[i][j+1]);
			if (i < q->blocks - 1) /* connect elements between blocks: */
				q->block[i][Q_BUFFER_SIZE - 1].next = &(q->block[i+1][0]);
		}
		q->block[q->blocks - 1][Q_BUFFER_SIZE - 1].next = NULL;
	}
	q->length = 0;
	q->head = NULL;
	return q;
}

void free_queue(QUEUE *q) {
	int i;
	if (q != NULL) {
		for (i = 0; i < q->blocks; i++)
			efree(q->block[i]); /* queue buffers */
		if (q->block != NULL)
			efree(q->block);
		efree(q);
	}
}

static Q_ELEMENT *get_free(QUEUE *q) {
	Q_ELEMENT *e;

	if (q->empty->next == NULL)
		enlarge_queue(q);
	e = q->empty;
	q->empty = q->empty->next;
	return e;
}

void enqueue(QUEUE *q, QUEUE_NODE *el, int n) {
/*
 * insert n elements in array el into the priority queue q
 */
	Q_ELEMENT *e, *where, *next;
	int i = 0, p;

	if (q == NULL || el == NULL || n <= 0)
		ErrMsg(ER_NULL, "enqueue");
	/* 
	 * first sort array el
	 */
	qsort(el, (size_t) n, sizeof(QUEUE_NODE), 
			(int CDECL (*)(const void *,const void *)) q->cmp);

	/*
 	 * and then merge them with the priority queue q
	 */

	/*
	 * find p, the number of elements in el[] that are smaller than
	 * the first element in the queue, if any.
	 * (Yes, we expect at least some of them to be closer than q->head)
	 */
	for (p = n; q->head != NULL && p > 0; p--)
		if (q->cmp(&(el[p-1]), &(q->head->el)) <= 0)
			break; /* out of this for loop */

	/* 
	 * put these p elements in order at the queue head: 
	 */
	for (i = p; i > 0; i--) {
		e = get_free(q);
		e->el = el[i-1];
		e->next = q->head;
		q->head = e;
	}
	q->length += p;
	n -= p;

	/*
	 * We might be done by now (n zero).
	 * Process the remaining elements:
	 */
	where = q->head; /* starting point for insertion: */
	next = where->next; /* the next element */
	el += p; /* start of the remaining elements */
	for (i = 0; i < n; i++) {
		e = get_free(q); /* get a free queue element */
		e->el = el[i]; /* copy contents: */
		/*
		 * unless where is the last element, shift a position in the queue
		 * when the element is larger than the insertion element e:
		 */
		while (next != NULL && q->cmp(&(e->el), &(next->el)) > 0) {
			where = next;
			next = where->next;
		}
		/*
		 * now next points either to NULL or to the first element smaller
		 * than or equal to e. So, insert e after where and before next:
   	 	 */
		e->next = next;
		where->next = e;
		/*
		 * the next element in el[] will should follow after e,
		 * so shift where one position to start looking after e:
		 */
		where = e;
	}
	q->length += n;
   	return;
}

QUEUE_NODE dequeue(QUEUE *q) {
	Q_ELEMENT *e;

	if (q->length == 0)
		ErrMsg(ER_NULL, "cannot dequeue empty queue");
	e = q->head; /* get first queue element */
	q->head = q->head->next; /* reset first to next */
	e->next = q->empty; /* put the dequeued element in the empty queue */
	q->empty = e;
	q->length--;
	return e->el;
}
