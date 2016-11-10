#ifndef DIMENSIONIZE_HH_
#define DIMENSIONIZE_HH_

template<class BaseDataType>
class Dimensionize{
private:
	std::vector<std::vector<BaseDataType>> DataPerDepth;
	std::vector<std::vector<BaseDataType>> DataPool;
	size_t max_depth;

	void DoIt(int myDepth, std::vector<BaseDataType> current_data){

		for(size_t i=0; i<DataPerDepth[myDepth-1].size(); i++){
			BaseDataType myData = DataPerDepth[myDepth-1][i];
			current_data[myDepth-1] = myData;

			if(myDepth == max_depth)
				DataPool.push_back(current_data);
			else
				DoIt(myDepth+1, current_data);
		}
	}

public:
	Dimensionize(size_t max_depth_, std::vector<std::vector<BaseDataType>> DataPerDepth_){
		max_depth = max_depth_;
		DataPerDepth = DataPerDepth_;

		if(DataPerDepth.size() != max_depth)
			throw std::runtime_error("DataPerDepth size mismatches max_depth !");
	}

	std::vector<std::vector<BaseDataType>> DoDimensionize(){
		DataPool.clear();
		std::vector<BaseDataType> current_data(max_depth);
		DoIt(1, current_data);
		return DataPool;
	}
};


#endif
