#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <math.h>
#include <time.h>

#define MAX(x, y) (((x) > (y)) ? (x) : (y))
#define MIN(x, y) (((x) < (y)) ? (x) : (y))

/* definition of the box in the grid according to the documentation */
struct box {
	int32_t id, y, x, height, width, top_cnt,
			bot_cnt, left_cnt, right_cnt;
	int32_t *top_ids, *bot_ids, *left_ids, *right_ids;
	double dsv;
};

/* global variables */
int g_boxes_cnt, g_rows_cnt, g_cols_cnt;
double g_max_dsv, g_min_dsv;
double g_affect_rate, g_epsilon, *g_tmp_dsv;
struct box *g_boxes;

/* functions declaration */
void read_neighbors(int *ids, int cnt);
void read_box(struct box *box);
bool is_converged();
double get_waat(struct box current_box);
void commit_changes();

void commit_changes()
{
	int i;
	for (i = 0; i < g_boxes_cnt; ++i)
		g_boxes[i].dsv = g_tmp_dsv[i];

}

bool is_converged()
{
	double min_dsv, max_dsv;
	int i;
	struct box current_box;
	bool res;

	res = false;
	min_dsv = g_boxes[0].dsv;
	max_dsv = g_boxes[0].dsv;

	for (i = 1; i < g_boxes_cnt; ++i) {
		current_box = g_boxes[i];
		max_dsv = MAX(current_box.dsv, max_dsv);
		min_dsv = MIN(current_box.dsv, min_dsv);
	}

#ifdef DEBUG
printf("diff: %lf, epsilon: %lf\n", max_dsv - min_dsv, epsilon * max_dsv);
#endif 

	if ((max_dsv - min_dsv) < g_epsilon * max_dsv) {
		g_max_dsv = max_dsv;
		g_min_dsv = min_dsv;
		res = true;
	}

	return res;
}

int get_contact_dist(struct box current, struct box neighbor, bool is_horizontal)
{
	int contact, cur_max, cur_min,
				 neighbor_max, neighbor_min;
	contact = 0;

	cur_min = !is_horizontal ? current.x : current.y;
	cur_max = !is_horizontal ? current.x + current.width : current.y + current.height;

	neighbor_min = !is_horizontal ? neighbor.x : neighbor.y;
	neighbor_max = !is_horizontal ? neighbor.x + neighbor.width : neighbor.y +
				   neighbor.height;

	contact = MIN(cur_max, neighbor_max) - MAX(neighbor_min, cur_min);

	return contact;

}

double get_waat(struct box current_box)
{
	double waat, sum;
	int perimeter, i, id, contact_dist;
	struct box neighbor;

	sum = 0;
	perimeter = (current_box.height + current_box.width) * 2;

	if (current_box.top_cnt == 0){
		sum += current_box.dsv * current_box.width;
	} else {
		for (i = 0; i < current_box.top_cnt; i++) {
			id = current_box.top_ids[i];
			neighbor = g_boxes[id];
			contact_dist = get_contact_dist(current_box, neighbor, false);
			sum += neighbor.dsv * contact_dist;
		}
	}

	if (current_box.bot_cnt == 0){
		sum += current_box.dsv * current_box.width;
	} else {
		for (i = 0; i < current_box.bot_cnt; i++) {
			id = current_box.bot_ids[i];
			neighbor = g_boxes[id];
			contact_dist = get_contact_dist(current_box, neighbor, false);
			sum += neighbor.dsv * contact_dist;
		}
	}

	if (current_box.right_cnt == 0){
		sum += current_box.dsv * current_box.height;
	} else {
		for (i = 0; i < current_box.right_cnt; i++) {
			id = current_box.right_ids[i];
			neighbor = g_boxes[id];
			contact_dist = get_contact_dist(current_box, neighbor, true);
			sum += neighbor.dsv * contact_dist;
		}
	}

	if (current_box.left_cnt == 0){
		sum += current_box.dsv * current_box.height;
	} else {
		for (i = 0; i < current_box.left_cnt; i++) {
			id = current_box.left_ids[i];
			neighbor = g_boxes[id];
			contact_dist = get_contact_dist(current_box, neighbor, true);
			sum += neighbor.dsv * contact_dist;
		}
	}

	waat = sum / perimeter;

	return waat;
}

void read_neighbors(int *ids, int cnt)
{
	int i;
	if (cnt == 0)
		*ids = -1;
	else
		for (i = 0; i < cnt; ++i) scanf(" %d", ids + i);
	scanf("\n");
}

void read_box(struct box *box)
{
	scanf("%d\n", &(box->id));
	scanf("%d %d %d %d\n",
			&(box->y), &(box->x),
			&(box->height), &(box->width));


	scanf("%d", &(box->top_cnt));
	box->top_ids = (int32_t *)malloc(box->top_cnt * sizeof(int32_t));
	read_neighbors(box->top_ids, box->top_cnt);


	scanf("%d", &(box->bot_cnt));
	box->bot_ids = (int32_t *)malloc(box->bot_cnt * sizeof(int32_t));
	read_neighbors(box->bot_ids, box->bot_cnt);

	scanf("%d", &(box->left_cnt));
	box->left_ids = (int32_t *)malloc(box->left_cnt * sizeof(int32_t));
	read_neighbors(box->left_ids, box->left_cnt);

	scanf("%d", &(box->right_cnt));
	box->right_ids = (int32_t *)malloc(box->right_cnt * sizeof(int32_t));
	read_neighbors(box->right_ids, box->right_cnt);

	scanf("%lf\n", &(box->dsv));
}

int main(int argc, char *argv[])
{
	int i, iter_cnt;
	double current_dsv, waat;
	time_t st, end;
	clock_t st_clock, end_clock;

	sscanf(argv[1], "%lf", &g_affect_rate);
	sscanf(argv[2], "%lf", &g_epsilon);


	/* parse the general grid info from the standard input */
	scanf("%d %d %d\n", &g_boxes_cnt, &g_rows_cnt, &g_cols_cnt);
	// printf("%d %d %d\n", boxes_cnt, rows_cnt, cols_cnt);
	g_boxes = (struct box *)malloc(g_boxes_cnt * sizeof(struct box));
	g_tmp_dsv = (double *)malloc(g_boxes_cnt * sizeof(double));

	/* parse each box in the grid */
	for (i = 0; i < g_boxes_cnt; i++)
		read_box(g_boxes + i);

	/* main modeling */
	iter_cnt = 0;
	time(&st);
	st_clock = clock();
	while(!is_converged()) {
		for (i = 0; i < g_boxes_cnt; ++i) {
			current_dsv = g_boxes[i].dsv;
			waat = get_waat(g_boxes[i]);
			g_tmp_dsv[i] = (get_waat(g_boxes[i]) - current_dsv)
						   * g_affect_rate + current_dsv;
#ifdef DEBUG
printf("%d    %lf    %lf    %lf\n", boxes[i].id, current_dsv, waat, tmp_dsv[i]);
#endif
		}
		commit_changes();
		iter_cnt++;
	}
	end_clock = clock();
	time(&end);

	printf("***************************************************************\n"
		   "dissipation converged in %d iterations,\n"
		   "	with max DSV = %lf and min DSV = %lf\n"
		   "	affect rate = %lf;     epsilon = %lf\n"
		   "elapsed convergence loop time (clock): %ld\n"
		   "elapsed convergence loop time (time): %ld\n"
		   "***************************************************************\n",
		   iter_cnt, g_max_dsv, g_min_dsv, g_affect_rate, g_epsilon,
		   end_clock - st_clock, end - st);

	return 0;
}
