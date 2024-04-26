from django import forms
from sugarDB.models import GLMFile

class GLMFileForm(forms.ModelForm):
    class Meta:
        model = GLMFile
        #fields = ('zipcode', 'file')
        fields = ('file',)
        widgets = {
            'zipcode': forms.TextInput(attrs={'placeholder': "Input the Zip Code Here"})
        }
        labels = {
            'zipcode': "Input the Zip Code Here",
            'file': "Upload a GLM File"
        }

    def clean_file(self):
        file = self.cleaned_data['file']
        if not file.name.endswith('.glm'):
            raise forms.ValidationError("The file is not a GLM file. Please upload a file with the .glm extension.")
        return file