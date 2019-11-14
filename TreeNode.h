//
// Created by Roman Ellerbrock on 2019-04-22.
//

#ifndef MCTDH_TREENODE_H
#define MCTDH_TREENODE_H
#include <utility>
#include <vector>
#include <iostream>

typedef vector<size_t> Path;

class NodeContent {
public:
    NodeContent() : f(0), n(0) {}
	explicit NodeContent(size_t n_) : f(0), n(n_) {}
    ~NodeContent()= default;

    void print(size_t indent = 0, std::ostream& os = std::cout)const {
        for (size_t i = 0; i < indent; ++i) {
            cout << "\t";
        }
        os << f << " " << n << std::endl;
    }

    size_t f, n;
};

class TreeNode {
public:

    /**
     * A simple tree-node class.
     */

    /// Rule of five-section (constructors & destructors)
    TreeNode() :up(nullptr), nextnodenr(0) {}
    explicit TreeNode(size_t n):up(nullptr), nextnodenr(0), content(n) {}
    ~TreeNode() = default;

    TreeNode(const TreeNode& node)
            : up(node.up), content(node.content), nextnodenr(node.nextnodenr){
        for (const TreeNode* child : node.down) {
            down.push_back(new TreeNode(*child));
        }
        for(TreeNode* child : down) {
            child->up = this;
        }
    }

    TreeNode(TreeNode&& node)noexcept
            : up(node.up), down(move(node.down)), content(node.content), nextnodenr(node.nextnodenr) {
        for (TreeNode* child : down) {
            child->up = this;
        }
    }

    TreeNode& operator=(const TreeNode& old) {
        TreeNode node(old);
        *this = std::move(node);
        return *this;
    }

    TreeNode& operator=(TreeNode&& old)noexcept {
        if (this == &old) {
            return *this;
        }

        up = old.up;
        down = move(old.down);
        content = old.content;
        nextnodenr = old.nextnodenr;

        for (TreeNode* child : down) {
            child->up = this;
        }
        return  *this;
    }

    /// Member functions
    TreeNode* NextNode() {
        TreeNode *result;
        if (nextnodenr < down.size()) {
            result = down[nextnodenr]->NextNode();
            if (result == down[nextnodenr]) {
                ++nextnodenr;
            }
        } else {
            nextnodenr = 0;
            result = this;
        }
        return result;
    }

    void push_back(const TreeNode& node) {
        down.emplace_back(new TreeNode(node));
    }

    size_t size()const { return down.size(); }

    const TreeNode& Down(size_t k)const {
        assert(k < down.size());
        return *down[k];
    }

    TreeNode& Down(size_t k) {
        assert(k < down.size());
        return *down[k];
    }

    void print(std::ostream& os = std::cout)const {
        content.print(Layer(), os);
        for (const TreeNode* node : down) {
            node->print(os);
        }
    }

    void GenInput(ostream& os = cout)const {
        // Indentation
        for (size_t i = 0; i < Layer() - 1; ++i) { os << "\t"; }
        if (IsLeaf()) {
            os << content.n << "\t" << "6" << "\t" << content.f << "\n";
        } else {
            os << content.n << "\t-" << size() << "\n";
        }
        for (TreeNode* child : down) {
            child->GenInput(os);
        }
    }

    bool IsLeaf()const { return (down.empty()); }
    bool IsRoot()const { return (up == nullptr); }

    void SetPath(const Path& newpath) {
        path = newpath;
        for (size_t k = 0; k < down.size(); ++k) {
            Path cpath = path;
            cpath.emplace_back(k);
            TreeNode& child = Down(k);
            child.SetPath(cpath);
        }
    }

	void SetPhysicalNodes(size_t& next) {
		if (IsLeaf()) {
			content.f = next;
			next++;
		}else {
			for (TreeNode* child : down) {
				child->SetPhysicalNodes(next);
			}
		}
	}

    void SetPhysicalNodesScatter(size_t& nextl, size_t& nextr, bool& left) {
        if (IsLeaf()) {
        	if (left) {
				content.f = nextl;
				nextl++;
        	} else {
				content.f = nextr;
				nextr++;
        	}
			left = !left;
        }else {
            for (TreeNode* child : down) {
                child->SetPhysicalNodesScatter(nextl, nextr, left);
            }
        }
    }

    void MakeRoot() {
        SetPath({0});
        up = nullptr;
    }

    size_t Layer()const { return path.size(); }

    NodeContent content;
private:
    TreeNode* up;
    std::vector<TreeNode*> down;
    int nextnodenr;
    Path path;
};


#endif //MCTDH_TREENODE_H
