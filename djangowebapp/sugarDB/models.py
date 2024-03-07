from django.db import models


class GLMFile(models.Model):
    zipcode = models.CharField(max_length=200, default='15213')
    file = models.FileField(upload_to='uploads/')
    uploaded_at = models.DateTimeField(auto_now_add=True)
