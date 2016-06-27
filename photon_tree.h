//
// Created by dc on 4/12/16.
//

#ifndef RAYTR_KDTREE_H
#define RAYTR_KDTREE_H

#include "common.h"

template <typename Info>
class PhotonTree
{
	typedef std::pair<Vec3, Info> Node;
	typedef std::vector<Node> NodeContainer;
	typedef typename NodeContainer::iterator NodeIterator;

	NodeContainer *container_;
	std::vector<double> lmax_k_, rmin_k_;

	void Build (NodeIterator s, NodeIterator e, int dim);
	template <typename Callback>
	void Query_ (NodeIterator s, NodeIterator e, const Vec3 &pos, double radius, Callback& callback, int dim);

	void DisKthNeighbor (NodeIterator s, NodeIterator e, const Vec3 &P, vector<double> &nbrs, int K, int dim);

public:

	int __counter;

	PhotonTree (NodeContainer *container);
	~PhotonTree ();
	template <typename Callback>
	void Query (const Vec3 &pos, double radius, Callback& callback);
	double DisKthNeighbor (const Vec3 &P, int K);
};

template <typename Info>
PhotonTree<Info>::PhotonTree (NodeContainer *container):
	container_(container)
{
	size_t size = container_->size();
	lmax_k_.resize(size);
	rmin_k_.resize(size);
	Build(container_->begin(), container_->end(), 0);
}
template <typename Info>
PhotonTree<Info>::~PhotonTree ()
{
	// delete container_;
}

template <typename Info>
template <typename Callback>
void PhotonTree<Info>::Query (const Vec3 &pos, double radius, Callback &callback)
{
	__counter = 0;
	Query_(container_->begin(), container_->end(), pos, radius, callback, 0);
}

template <typename Info>
template <typename Callback>
void PhotonTree<Info>::Query_ (PhotonTree<Info>::NodeIterator s, PhotonTree<Info>::NodeIterator e,
						   const Vec3 &pos, double radius, Callback &callback, int dim)
{
	if (s >= e)
		return;

	auto m = s + ((e - s) >> 1);
	if ((m->first - pos).sqnorm() <= sqr(radius))
	{
		callback(m->first, m->second);
		++__counter;
	}

	if (s + 1 < e)
	{
		int offs = m - container_->begin();
		if (pos[dim] - radius <= lmax_k_[offs])
			Query_(s, m, pos, radius, callback, (dim + 1) % 3);

		if (pos[dim] + radius >= rmin_k_[offs])
			Query_(m + 1, e, pos, radius, callback, (dim + 1) % 3);
	}
}

template <class Info>
double PhotonTree<Info>::DisKthNeighbor (const Vec3 &P, int K)
{
	vector<double> candidates(K, inf);
	DisKthNeighbor(container_->begin(), container_->end(), P, candidates, K + 1, 0);
	for (auto it = candidates.rbegin(); it != candidates.rend(); ++it)
	{
		if (*it < inf * 0.5)
			return sqrt(*it);
	}
	assert(false);
}

template <class Info>
void PhotonTree<Info>::DisKthNeighbor (NodeIterator s, NodeIterator e, const Vec3 &P, vector<double> &nbrs, int K,
									   int dim)
{
	if (s >= e)
		return;

	if (s + 1 == e)
	{
		auto sqdis = (s->first - P).sqnorm();
		if (sqdis < nbrs.back())
			*nbrs.rbegin() = sqdis;
	}
	else
	{
		auto m = s + ((e - s) >> 1);
		int offs = m - container_->begin();
		if (P[dim] <= lmax_k_[offs])
		{
			DisKthNeighbor(s, m, P, nbrs, K, (dim + 1) % 3);
			if (sqr(P[dim] - rmin_k_[offs]) < nbrs.back())
				DisKthNeighbor(m, e, P, nbrs, K, (dim + 1) % 3);
		}
		else
		{
			DisKthNeighbor(m, e, P, nbrs, K, (dim + 1) % 3);
			if (sqr(P[dim] - lmax_k_[offs]) < nbrs.back())
				DisKthNeighbor(s, m, P, nbrs, K, (dim + 1) % 3);
		}
	}
}

template <typename Info>
void PhotonTree<Info>::Build (PhotonTree<Info>::NodeIterator s, PhotonTree<Info>::NodeIterator e, int dim)
{
	if (s >= e)
		return;

	auto mms = (e - s) >> 1;
	auto dcmp = [dim](const Node &a, const Node &b){ return a.first[dim] + eps < b.first[dim]; };
	GP::nth_element(s, s + mms, e, dcmp);
	auto m = s + mms;
	ptrdiff_t offset = m - container_->begin();
	lmax_k_[offset] = m == s ? m->first[dim] : GP::max_element(s, m, dcmp)->first[dim];
	rmin_k_[offset] = m + 1 == e ? m->first[dim] : GP::min_element(m + 1, e, dcmp)->first[dim];
	Build(s, m, (dim + 1) % 3);
	Build(m + 1, e, (dim + 1) % 3);
}


#endif //RAYTR_KDTREE_H
