// small_vector.h
#ifndef SMALL_VECTOR_H
#define SMALL_VECTOR_H

#include <vector>
#include <stdexcept>
#include <cstddef>
#include <cassert>

template <typename T>
class SmallVector {
private:
	static const size_t BLOCK_SIZE = 8192; // 每个块的容量，可以根据需要调整
	std::vector<T*> blocks;                  // 存储指向每个块的指针
	size_t vec_size;                         // 当前元素数量

	// 分配一个新的块
	void allocate_block() {
		T* new_block = new T[BLOCK_SIZE];
		blocks.push_back(new_block);
	}

public:
	// 构造函数
	SmallVector() : vec_size(0) {}

	// 析构函数
	~SmallVector() {
		for (auto block : blocks) {
			delete[] block;
		}
	}

	// 禁用拷贝构造函数和拷贝赋值运算符
	SmallVector(const SmallVector&) = delete;
	SmallVector& operator=(const SmallVector&) = delete;

	// 启用移动构造函数和移动赋值运算符
	SmallVector(SmallVector&& other) noexcept
		: blocks(std::move(other.blocks)), vec_size(other.vec_size) {
		other.vec_size = 0;
	}

	SmallVector& operator=(SmallVector&& other) noexcept {
		if (this != &other) {
			for (auto block : blocks) {
				delete[] block;
			}
			blocks = std::move(other.blocks);
			vec_size = other.vec_size;
			other.vec_size = 0;
		}
		return *this;
	}

	// 添加元素到容器末尾
	void push_back(const T& value) {
		size_t block_index = vec_size / BLOCK_SIZE;
		size_t index_in_block = vec_size % BLOCK_SIZE;

		// 如果当前块不存在或已满，分配一个新的块
		if (block_index >= blocks.size()) {
			allocate_block();
		}

		// 在当前块中添加元素
		blocks[block_index][index_in_block] = value;
		++vec_size;
	}

	// 重置容器大小，并根据需要填充新元素
	void resize(size_t new_size, const T& value = T()) {
		if (new_size > vec_size) {
			// 确保有足够的块
			size_t required_blocks = (new_size + BLOCK_SIZE - 1) / BLOCK_SIZE;
			while (blocks.size() < required_blocks) {
				allocate_block();
			}
			// 填充新元素
			for (size_t i = vec_size; i < new_size; ++i) {
				size_t block_index = i / BLOCK_SIZE;
				size_t index_in_block = i % BLOCK_SIZE;
				blocks[block_index][index_in_block] = value;
			}
		}
		// 如果缩小大小，只需调整 vec_size
		vec_size = new_size;
	}

	// 预留容量
	void reserve(size_t new_capacity) {
		if (new_capacity > vec_size) {
			size_t required_blocks = (new_capacity + BLOCK_SIZE - 1) / BLOCK_SIZE;
			while (blocks.size() < required_blocks) {
				allocate_block();
			}
		}
	}

	// 获取当前大小
	size_t size() const {
		return vec_size;
	}

	// 获取当前容量
	size_t get_capacity() const {
		return blocks.size() * BLOCK_SIZE;
	}

	// 重载下标操作符（非 const 版本）
	T& operator[](size_t index) {
		if (index >= vec_size) {
			throw std::out_of_range("Index out of range");
		}
		size_t block_index = index / BLOCK_SIZE;
		size_t index_in_block = index % BLOCK_SIZE;
		return blocks[block_index][index_in_block];
	}

	// 重载下标操作符（const 版本）
	const T& operator[](size_t index) const {
		if (index >= vec_size) {
			throw std::out_of_range("Index out of range");
		}
		size_t block_index = index / BLOCK_SIZE;
		size_t index_in_block = index % BLOCK_SIZE;
		return blocks[block_index][index_in_block];
	}
};

#endif // SMALL_VECTOR_H
