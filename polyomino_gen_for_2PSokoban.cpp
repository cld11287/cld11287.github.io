/* 再帰を使わないポリオミノ生成 + 双庫番に適したものへの絞り込み */
/* 最終更新日: 2023/10/4 */
/* ハードコードする変数：AREA_SIZE */
/* mathmasterzachさんがdiscordで送ってくれた再帰あり・地形検出なしのコードの使い回しです */
/* http://www.mathmasterzach.com/ */
/* neighborsを右上下左の順にする */

#include <array>
#include <algorithm>//sort
#include <bitset>
#include <fstream>
#include <iostream>
#include <string>
#include <time.h>
#include <tuple>
#include <vector>

/*
面積 1 2 3 4  5  6   7   8    9    10    11     12     13      14      15       16
両面 1 1 2 5  12 35  108 369  1285 4655  17073  63600  238591  901971  3426576
有向 1 2 6 19 63 216 760 2725 9910 36446 135268 505861 1903890 7204874 27394666
*/

constexpr int AREA_SIZE = 17;
std::string filename = "../poly_" + std::to_string(AREA_SIZE) + ".txt";

std::ofstream afile;

struct Square { int x, y; };
/* zachはSquareより２次元配列のが速いって言ってたけどそのインデックスを格納したらSquareやん */
/* エイゴワカリマセン */
typedef std::array<Square, AREA_SIZE> Pomino;

inline bool operator ==(const Square& lhs, const Square& rhs) {
	return lhs.x == rhs.x && lhs.y == rhs.y;
}

inline bool operator !=(const Square& lhs, const Square& rhs) {
	return lhs.x != rhs.x || lhs.y != rhs.y;
}

inline bool operator <(const Square& lhs, const Square& rhs) {
	return lhs.x < rhs.x || (lhs.x == rhs.x && lhs.y < rhs.y);
}

/* 回転したり反転したりして重複を確認する際に座標を合わせるためのやつ */
Pomino clean_ordering(const Pomino& vec) {
	Square min_sq = { 0, 0 };
	for (const Square& s : vec) min_sq = { std::min(min_sq.x, s.x), std::min(min_sq.y, s.y) };
	Pomino tmp{};
	for (int i = 0; i < AREA_SIZE; ++i) tmp[i] = { vec[i].x - min_sq.x, vec[i].y - min_sq.y };
	return tmp;
}

/* ファイルに書き込む用 */
std::string print_polyomino(Pomino& vec) {
	Square min_sq = { 0, 0 }, max_sq = { 0, 0 };
	for (const Square& s : vec) {
		min_sq = { std::min(min_sq.x, s.x), std::min(min_sq.y, s.y) };
		max_sq = { std::max(max_sq.x, s.x), std::max(max_sq.y, s.y) };
	}
	std::string str{ "" };
	for (int y = max_sq.y; y >= min_sq.y; --y) {
		for (int x = min_sq.x; x <= max_sq.x; ++x) {
			bool found = false;
			for (const Square& s : vec) if (s.x == x && s.y == y) found = true;
			if (found) str += '-';
			else str += '#';
		}
		// replace this "\n" with "|" to keep the whole polyomino on one line
		if (y != min_sq.y) str += '|';
	}
	str += '\n';
	return str;
}

/* 出力用 */
void printout_polyomino(const Pomino& vec) {
	Square min_sq = { 0, 0 }, max_sq = { 0, 0 };
	for (const Square& s : vec) {
		min_sq = { std::min(min_sq.x, s.x), std::min(min_sq.y, s.y) };
		max_sq = { std::max(max_sq.x, s.x), std::max(max_sq.y, s.y) };
	}
	std::string str{ "" };
	for (int y = max_sq.y; y >= min_sq.y; --y) {
		for (int x = min_sq.x; x <= max_sq.x; ++x) {
			bool found = false;
			for (const Square& s : vec) if (s.x == x && s.y == y) found = true;
			if (found) str += '-';
			else str += '#';
		}
		// replace this "\n" with "|" to keep the whole polyomino on one line
		if (y != min_sq.y) str += '\n';
	}
	str += "\n\n";
	std::cout << str;
}

/*
周囲３マスを無効にする　HH = AREA_SIZE / 2
AREA_SIZE = 6 or 7 のとき
###########
###########
###########
###-----###
###-----###
#####o--###
###########
###########
###########
*/
/* 検出する地形のセルの基準点からの距離の最大値に応じて変化させる */
constexpr int WW = (AREA_SIZE / 2) * 2 - 1 + 6;
constexpr int bset_size = WW * (AREA_SIZE / 2 + 6);
constexpr int ORIGIN = 3 + (AREA_SIZE / 2) + WW * 3 - 1;

/* 近隣のセル　インデックスを３から引けば時計回りに90度回せる */
std::array<int, 4> nearby { 1, -1, WW, -WW };
std::bitset<bset_size> boundout_master = 0, boundout = 0;
std::vector<std::tuple<std::bitset<bset_size>, std::bitset<bset_size>, int, int>> detector_last{};
std::vector<std::array<std::bitset<bset_size>, 2>> detector_middle{};

void initiate() {
	for (int a = 0; a < ORIGIN; ++a) boundout_master.set(a);
	for (int a = 3; a < 3 + AREA_SIZE / 2; ++a) {
		for (int b = 0; b < 3; ++b) boundout_master.set(a * WW + b);
		for (int b = WW - 3; b < WW; ++b) boundout_master.set(a * WW + b);
	}
	for (int a = 3 + AREA_SIZE / 2; a < 6 + AREA_SIZE / 2; ++a) {
		for (int b = 0; b < WW; ++b) boundout_master.set(a * WW + b);
	}
	/*
	704  -  #-- #--
	3 1 #-# --# ---
	625  -      --#
	*/
	detector_last.resize(8);
	std::vector<int> dir{ WW, 1, -WW, -1, WW + 1, -WW + 1, -WW - 1, WW - 1 };
	std::vector<std::vector<int>> index_vec{};
	index_vec.resize(16);
	/* pattern 1 */
	index_vec[0] = { 0, 2 };
	index_vec[8] = { 1, 3 };
	index_vec[1] = { 1, 3 };
	index_vec[9] = { 0, 2 };
	/* 2　1と区別するために隣にstartが来ないようにする */
	index_vec[2] = { 4, 3, 0 };
	index_vec[10] = { 7, 1 };
	index_vec[3] = { 7, 1, 0 };
	index_vec[11] = { 4, 3 };
	index_vec[4] = { 4, 2, 1 };
	index_vec[12] = { 5, 0 };
	index_vec[5] = { 5, 0, 1 };
	index_vec[13] = { 4, 2 };
	/* 3 */
	index_vec[6] = { 4, 6, 0, 1, 2, 3 };
	index_vec[14] = { 7, 5 };
	index_vec[7] = { 7, 5, 0, 1, 2, 3 };
	index_vec[15] = { 4, 6 };
	for (int a = 0; a < 8; ++a) {
		for (const int& b : index_vec[a]) std::get<0>(detector_last[a]).set(ORIGIN + dir[b]);
		for (const int& b : index_vec[a + 8]) std::get<1>(detector_last[a]).set(ORIGIN + dir[b]);
		std::get<2>(detector_last[a]) = ORIGIN + dir[index_vec[a][0]];
		std::get<3>(detector_last[a]) = ORIGIN + dir[index_vec[a][1]];
	}
	/* ここからmiddle */
	std::array<std::array<std::vector<Square>, 2>, 38> dir_table{};
	/* 蓑のセルを基準にする　７個目のパターンは角に２マス付いた地形を排除してしまうけど別にいいや */
	/* #   #    #     # #  ##    #     #      #
		  #@-# #@ -# #@        #  -#   -# # @ -#
	  #@#   #     #   ### #@-# #@  # #@    ##
	   #                   ##    #    #          */
	dir_table[0][0] = {};
	dir_table[0][1] = { { -1, 0 }, { 0, 2 }, { 0, -1 }, { 1, 0 } };
	dir_table[4][0] = { { 1, 0 } };
	dir_table[4][1] = { { -1, 0 }, { 0, 1 }, { 1, -1 }, { 2, 0 } };
	dir_table[8][0] = { { 2, 0 } };
	dir_table[8][1] = { { -1, 0 }, { 0, 1 }, { 2, -1 }, { 3, 0 } };
	dir_table[12][0] = {};
	dir_table[12][1] = { { -1, 0 }, { 0, 1 }, { 0, -1 }, { 1, -1 }, { 2, 1 }, { 2, -1 } };
	dir_table[20][0] = { { 1, 0 } };
	dir_table[20][1] = { { -1, 0 }, { 0, 2 }, { 0, -1 }, { 1, 2 }, { 1, -1 }, { 2, 0 } };
	dir_table[24][0] = { { 2, 1 } };
	dir_table[24][1] = { { -1, 0 }, { -1, 1 }, { 1, -1 }, { 1, 2 }, { 3, 0 }, { 3, 1 } };
	dir_table[28][0] = { { 1, 1 } };
	dir_table[28][1] = { { -1, 0 }, { 0, -1 }, { 1, 2 }, { 2, 1 } };
	dir_table[30][0] = { { 2, 0 } };
	dir_table[30][1] = { { -2, 0 }, { -1, -1 }, { 0, -1 }, { 2, 1 }, { 3, 0 } };
	for (int a : { 12, 30 }) {
		for (int b = 0; b < 2; ++b) {
			for (int i = 0; i < 7; ++i) {
				/* 水平方向に反転 */
				if (i == 3) for (const Square& freq : dir_table[a + i][b]) dir_table[a + i + 1][b].push_back({ -freq.x, freq.y });
				/* 時計回りに90°回転 */
				else for (const Square& freq : dir_table[a + i][b]) dir_table[a + i + 1][b].push_back({ freq.y, -freq.x });
			}
		}
	}
	/* 線対称の場合 */
	for (int a : { 0, 20 }) {
		for (int b = 0; b < 2; ++b) {
			for (int i = 0; i < 3; ++i) {
				/* 時計回りに90°回転 */
				for (const Square& freq : dir_table[a + i][b]) dir_table[a + i + 1][b].push_back({ freq.y, -freq.x });
			}
		}
	}
	/* 点対称の場合１ */
	for (int a : { 4, 8, 24 }) {
		for (int b = 0; b < 2; ++b) {
			for (int i = 0; i < 3; ++i) {
				/* 時計回りに90°回転 */
				if (i == 1) for (const Square& freq : dir_table[a + i][b]) dir_table[a + i + 1][b].push_back({ freq.y, -freq.x });
				/* 水平方向に反転 */
				else for (const Square& freq : dir_table[a + i][b]) dir_table[a + i + 1][b].push_back({ -freq.x, freq.y });
			}
		}
	}
	/* 点対称の場合２ */
	for (int a : { 28 }) {
		for (int b = 0; b < 2; ++b) {
			/* 時計回りに90°回転 */
			for (const Square& freq : dir_table[a][b]) dir_table[a + 1][b].push_back({ freq.y, -freq.x });
		}
	}
	detector_middle.resize(dir_table.size());
	for (int a = 0; a < dir_table.size(); ++a) for (int b = 0; b < 2; ++b) {
		for (const Square& sq : dir_table[a][b]) {
			detector_middle[a][b].set(ORIGIN + sq.x + WW * sq.y);
		}
	}
}

bool is_deficient_middle(const Pomino& board, const std::vector<Square>& disabled, const int& poly_size) {
	boundout = boundout_master;
	for (const Square& a : disabled) boundout.set(ORIGIN + a.x + WW * a.y);
	std::bitset<bset_size> board_bs = 0;
	for (const Square& a : board) board_bs.set(ORIGIN + a.x + WW * a.y);
	int rel = 0;
	for (int a = 0; a < poly_size; ++a) {
		rel = board[a].x + WW * board[a].y;
		for (const std::array<std::bitset<bset_size>, 2>&b : detector_middle) {
			/* intを返す場合、無効セルだけの地形を検出した場合はそのインデックスを返さなければいけない（非常にめんどくさい） */
			/* lastは地形が排反だがmiddleは排反ではないのでbreakできない　少しいじったらいけそうだけど */
			if ((~board_bs & (b[0] << rel)).none() && (~boundout & (b[1] << rel)).none()) return true;
		}
		int m = 3, n = 0;
		for (const int& b : nearby) {
			--m;
			if (!boundout[ORIGIN + rel + b]) {
				++n;
				m = 3;
			}
		}
		if (n != 1) continue;
		int point = ORIGIN + rel;
		while (true) {
			if (!boundout[point + nearby[3 - m]] && !boundout[point - nearby[3 - m]]) break;
			if (boundout[point + nearby[m]]) return true;
			else if (!board_bs[point + nearby[m]]) break;
			point += nearby[m];
			if (point < WW * 4) return true;
		}
	}
	return false;
}

int is_separated(const int& start, const int& end, const int& marked, const std::bitset<bset_size>& board_bs) {
	bool flag = false;
	int flood_size = 1;
	std::bitset<bset_size> frontier = 0, new_frontier = 0, searched = 0;
	new_frontier.set(start);
	searched.set(start);
	searched.set(marked);
	while (new_frontier.any()) {
		frontier = new_frontier;
		new_frontier = 0;
		/* frontierを走査 */
		for (int a = ORIGIN; a < bset_size - WW * 3 - 3; ++a) {
			if (!frontier[a]) continue;
			for (const int& b : nearby) {
				if ((a == start) && (a + b == marked)) flag = true;
				if (!board_bs[a + b] || searched[a + b]) continue;
				if (a + b == end) return 0;
				++flood_size;
				new_frontier.set(a + b);
			}
			/* あるマスから見た４方向は重ならないのでここで */
			searched |= new_frontier;
		}
	}
	/* 長さ２の袋小路の場合（パターン１）は除きたい　途中で検出する場合はAREA_SIZEではなくboard_bs.count()にする */
	if (flag && (flood_size == 1 || flood_size == board_bs.count() - 2)) return 0;
	return flood_size;
}

bool is_deficient_last(const Pomino& board) {
	std::bitset<bset_size> board_bs = 0;
	for (const Square& a : board) board_bs.set(ORIGIN + a.x + WW * a.y);
	int rel = 0;
	for (int a = 0; a < AREA_SIZE; ++a) {
		rel = board[a].x + WW * board[a].y;
		for (const std::array<std::bitset<bset_size>, 2>&b : detector_middle) {
			if ((~board_bs & (b[0] << rel)).none() && (board_bs & (b[1] << rel)).none()) return true;
		}
		int m = 3, n = 0;
		for (const int& b : nearby) {
			--m;
			if (board_bs[ORIGIN + rel + b]) {
				++n;
				m = 3;
			}
		}
		if (n == 1) {
			int point = ORIGIN + rel;
			while (true) {
				if (board_bs[point + nearby[3 - m]] && board_bs[point - nearby[3 - m]]) break;
				if (!board_bs[point + nearby[m]]) return true;
				point += nearby[m];
			}
		}
		for (const std::tuple<std::bitset<bset_size>, std::bitset<bset_size>, int, int>& b : detector_last) {
			if ((~board_bs & (std::get<0>(b) << rel)).none() && (board_bs & (std::get<1>(b) << rel)).none()) {
				if (is_separated(std::get<2>(b) + rel, std::get<3>(b) + rel, ORIGIN + rel, board_bs) > 0) return true;
				break;
			}
		}
	}
	return false;
}

int suitable_poly_count = 0;

void finalize_polyomino(const Pomino& vec) {
	Pomino orig = clean_ordering(vec), test = orig, tmp{};
	std::sort(orig.begin(), orig.end());
	/* もしbitsetを使うなら回転するのではなく隅の８方向からの辞書順を確かめればいい */
	for (int i = 0; i < 7; ++i) {
		if (i == 3) {
			/* 水平方向に反転 */
			for (int j = 0; j < AREA_SIZE; ++j) tmp[j] = { -test[j].x, test[j].y };
			test = clean_ordering(tmp);
		}
		else {
			/* 時計回りに90°回転 */
			for (int j = 0; j < AREA_SIZE; ++j) tmp[j] = { test[j].y, -test[j].x };
			test = clean_ordering(tmp);
		}
		std::sort(test.begin(), test.end());
		/* 辞書的順序で絞り込む　回転してboundout_masterからはみ出る蓑が優先されてもよい */
		for (int i = 0; i < AREA_SIZE; ++i) if (test[i] != orig[i]) {
			if (test[i] < orig[i]) return;
			break;
		}
	}
	if (is_deficient_last(vec)) return;
	//printout_polyomino(orig);
	afile << print_polyomino(orig);
	++suitable_poly_count;
}

int generated_count = 0, detect_count = 0;

/* y <= -1 || (y == 0 && x <= -1) を無効セルにしておく */
/* 蓑の最新のセルiの周囲のセルにi+1から番号を振り、無効セルでない最も古いセルを追加する */
/* 蓑が最大になったら蓑の最も新しいセルを無効セルにする */
/* 蓑の周りが全て無効セルになったら蓑で最も新しいセルより新しい全ての無効セルを解放し、蓑の最も新しいセルを無効セルにする */
/* 解放する際、蓑で二番目に新しいセルの周囲のセルで蓑の最も新しいセルより新しいセルは番号ごと空ける */
/* 生成自体では無効セルを直接用いていない（cursorで解決している）が、地形検出のために結局使う */
void enumerate() {
	/* frontier_sizeは境界セルのサイズではなくfrontier_vecのサイズです　size()で計算してもほとんど速度変わらなかったけどこっちが好き */
	/* cursorはindexそのものなのでsizeより1小さい */
	int cursor = 0, poly_size = 0, frontier_size = 1;
	std::array<Square, 4> neighbors{};
	neighbors[0] = { 1, 0 };
	neighbors[1] = { -1, 0 };
	neighbors[2] = { 0, 1 };
	neighbors[3] = { 0, -1 };
	/* 何が起こってるかよく分からんけど動いてるのでヨシ！ */
	std::array<std::array<int, 2>, AREA_SIZE> cursor_vec{};
	Pomino poly_vec{};
	/* 無効セルを格納する　無効セルよりも新しいセルが蓑に含まれていることがあるので注意 */
	std::vector<Square> disabled{};
	/* frontier_vecって名前だけどポリオミノとその周辺の無効含む全てのセルが含まれる */
	/* arrayにしたら逆に遅くなった　push_backよりも値を参照する方が遅いらしい */
	std::vector<Square> frontier_vec{ { 0, 0 } };
	while (true) {
		if (cursor < frontier_size) {
			/* 追加　両辺のインデックスは当然常には一致しない */
			poly_vec[poly_size] = frontier_vec[cursor];
			++cursor;
			if (poly_size == AREA_SIZE - 1) {
				++generated_count;
				finalize_polyomino(poly_vec);
				disabled.push_back(poly_vec[poly_size]);
				poly_vec[poly_size] = {};
				continue;
			}
			/* cursorは次に追加する境界セルになっている */
			cursor_vec[poly_size] = { frontier_size, cursor };
			++poly_size;
			/* frontierを拡張 */
			for (const Square& a : neighbors) {
				Square tmp_sq = { a.x + poly_vec[poly_size - 1].x, a.y + poly_vec[poly_size - 1].y };
				if (!boundout_master[ORIGIN + tmp_sq.x + tmp_sq.y * WW] && std::find(frontier_vec.begin(), frontier_vec.end(), tmp_sq) == frontier_vec.end()) {
					frontier_vec.push_back(tmp_sq);
					++frontier_size;
				}
			}
		}
		else {
			if (poly_size == 1) return;
			for (int i = frontier_size - cursor_vec[poly_size - 1][1]; i > 0; --i) disabled.pop_back();
			/* 蓑から１つ減らした際の最後の境界セルの１つ上まで削除する */
			/* resizeを使いたい */
			frontier_vec.resize(cursor_vec[poly_size - 1][0]);
			frontier_size = cursor_vec[poly_size - 1][0];
			/* 下２行は２連続でelseに入ったときしか使わない */
			poly_vec[poly_size] = {};
			cursor_vec[poly_size] = {};
			--poly_size;
			cursor = cursor_vec[poly_size][1];
			disabled.push_back(poly_vec[poly_size]);
			/* 下の文が無いと蓑と無効セルがここで被るけど検出以外では問題ないらしい */
			poly_vec[poly_size] = {};
			++detect_count;
			/* 境界セルから地形が成立する最小のセルを割り出せたらfrontierを消す都度調べなくて済む気がしたけど無理そう */
			if (cursor < frontier_size && is_deficient_middle(poly_vec, disabled, poly_size)) {
				cursor = frontier_size;
				//if (cursor_vec[poly_size][1] - (int)disabled.size() != poly_size) std::cout << "!\n";
				disabled.resize(frontier_size - poly_size);
			}
		}
	}
}

int main() {
	std::cout << "Search for polyominoes of size: " << AREA_SIZE << "\n";
	std::ifstream bfile;
	bfile.open(filename);
	std::string tmp_st;
	if (bfile.is_open() && std::getline(bfile, tmp_st)) {
		std::cout << "Are you sure you want to overwrite the previous output file?\nYes : press Y\nNo : press N" << std::endl;
		char tmp;
		while (true) {
			std::cin >> tmp;
			if (tmp == 'y' || tmp == 'Y') break;
			if (tmp == 'n' || tmp == 'N') exit(1);
		}
	}
	else std::cout << "Created a new file : " << filename << std::endl;
	bfile.close();
	std::cout << "Calculation started" << std::endl;
	clock_t tstart = clock();
	initiate();
	afile.open(filename);
	enumerate();
	afile.close();
	std::cout << "suitable: " << suitable_poly_count << " out of: " << generated_count << " have inspected: " << detect_count;
	std::cout << "\n\n" << (double)(clock() - tstart) / CLOCKS_PER_SEC << "sec" << "\nend";
	return 0;
}
